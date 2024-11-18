

### Smart Pointers in C++: A Beginner's Guide for Mesh Processing

If you’ve worked with pointers in other languages, you’re probably used to managing memory manually—allocating memory for new objects and making sure to free it when you’re done. In C++, this has traditionally been done with raw pointers, but they come with a big downside: **memory management is entirely your responsibility**. If you forget to free memory or accidentally use freed memory, your program could crash, or you could introduce subtle bugs, like memory leaks.

Luckily, C++ provides **smart pointers**, which help manage memory for you. Smart pointers take care of cleaning up memory when it’s no longer needed, making your code safer and easier to work with. The three most common types are `std::unique_ptr`, `std::shared_ptr`, and `std::weak_ptr`. Let’s dive into what each of these does and how they differ from regular pointers.

---

#### 1. Raw Pointers: The Old School Way

Before smart pointers, C++ developers used **raw pointers** to manage memory. For example:

```cpp
Mesh* mesh = new Mesh();
```

This allocates memory for a `Mesh` object and returns a pointer to that memory. However, you must **manually delete** this object when you no longer need it:

```cpp
delete mesh;
```

If you forget to call `delete`, your program will have a **memory leak**—the memory that `mesh` points to won’t be freed. Also, if you try to use `mesh` after calling `delete`, you’ll be working with a **dangling pointer**, which can lead to undefined behavior.

---

#### 2. `std::unique_ptr`: Exclusive Ownership

`std::unique_ptr` is the simplest type of smart pointer. It enforces **exclusive ownership** of a resource. This means only one `std::unique_ptr` can point to a given object at a time, so it’s perfect when you know that an object should only ever be owned by one part of your program.

Here’s an example:

```cpp
std::unique_ptr<Mesh> mesh = std::make_unique<Mesh>();
```

Unlike raw pointers, you don’t have to call `delete` on this pointer. When `mesh` goes out of scope (e.g., when the function it’s in returns), C++ will automatically free the memory for you. If you try to copy a `std::unique_ptr`, you’ll get a compiler error, because it’s designed to prevent multiple owners.

```cpp
std::unique_ptr<Mesh> mesh2 = mesh;  // ERROR: can't copy a unique_ptr
```

However, you can **transfer ownership** using `std::move()`:

```cpp
std::unique_ptr<Mesh> mesh2 = std::move(mesh);  // mesh is now empty
```
In our code, it is for example convenient to use `std::unique_ptr` to reference an instance of a mesh class. 

---

#### 3. `std::shared_ptr`: Shared Ownership

Sometimes, you need multiple parts of your program to share ownership of a resource. In mesh processing, for instance, you might have multiple algorithms operating on the same mesh. For this, C++ provides `std::shared_ptr`.

Here’s an example:

```cpp
std::shared_ptr<Mesh> mesh1 = std::make_shared<Mesh>();
std::shared_ptr<Mesh> mesh2 = mesh1;  // both mesh1 and mesh2 own the object
```

With `std::shared_ptr`, an object’s memory is freed only when the **last** `shared_ptr` that points to it is destroyed. C++ keeps track of how many `std::shared_ptr`s refer to the object with something called a **reference count**. When the count drops to zero, the memory is released.

However, this convenience comes with some cost: `std::shared_ptr` is a little slower than `std::unique_ptr`, since it has to maintain and check the reference count.

---

#### 4. `std::weak_ptr`: Non-owning Reference

There are cases when you need a pointer to an object without affecting its reference count. Enter `std::weak_ptr`. It’s like a non-owning observer of an object that is managed by `std::shared_ptr`.

Imagine you have two objects that reference each other, such as a `Mesh` and an `Edge`. If they both use `std::shared_ptr`, you can end up with a **circular reference**, where each object holds onto the other indefinitely, preventing their memory from being freed. This is where `std::weak_ptr` shines.

```cpp
std::shared_ptr<HalfEdge> hedge = std::make_shared<HalfEdge>();
std::weak_ptr<HalfEdge> weakhedge = hedge;  // does not increase the reference count
```

`std::weak_ptr` doesn’t manage the memory directly, so it won’t keep the object alive. If all `std::shared_ptr`s to the object are destroyed, the object’s memory will be freed, even if a `std::weak_ptr` still exists.

When you want to access the object, you need to **lock** the `std::weak_ptr` to create a temporary `std::shared_ptr`.
In other words `std::weak_ptr` is a smart pointer that holds a non-owning reference to
an object managed by `std::shared_ptr`. 

This is for example what we are doing in the `MeshParts.h` when we want to set the pointers for flipped half edges and next half edges. We basically check if the object that the pointer points to exists, if yes, the corresponding shared pointer is returned, otherwise we throw a runtime error:

```cpp
struct HalfEdge{ 
    [...]

    HalfEdgePtr getFlipHalfEdge() const
    {
        if (flip.expired())
        {
            throw std::runtime_error("Flip Half Edge is deleted");
        }
        else
        {
            return flip.lock();
        }
    }

    HalfEdgePtr getNextHalfEdge() const
    {

        if (next.expired())
        {
            throw std::runtime_error("Next Half Edge is deleted");
        }
        else
        {
            return next.lock();
        }
    }

    PrimalFacePtr getPrimalFace() const
    {
        if (primal_face.expired())
        {
            throw std::runtime_error("Primal Face is deleted");
        }
        else
        {
            return primal_face.lock();
        }
    }
    [...]
}
```
For instance, if we would use `std::shared_ptr` to define the half-edge datastructure, we would run into an issue due to cyclic references..

---

#### Why Use Smart Pointers?

Smart pointers make your C++ code safer and easier to work with. They help you avoid:

- **Memory leaks**: Smart pointers automatically free memory when it’s no longer needed.
- **Dangling pointers**: Smart pointers become null after their object is destroyed.
- **Circular references**: Using `std::weak_ptr` with `std::shared_ptr` avoids keeping objects alive unnecessarily.

By using smart pointers, you don’t have to worry as much about manually managing memory, which reduces bugs and improves your code’s stability. In mesh processing, where data structures can get complex and interdependent, using smart pointers will simplify your life while ensuring that resources are properly managed.




### Notes on Accessing, Referencing, and Dereferencing

- **Referencing**: 
  - You create a smart pointer like this: 
    ```cpp
    std::unique_ptr<Mesh> meshPtr = std::make_unique<Mesh>();
    ```
  - Smart pointers also support normal pointer-like operations such as referencing objects (`&object`) to get their memory address, but generally, you manage the pointer through the smart pointer itself.

- **Dereferencing**:
  - To access the object stored in a smart pointer, you can dereference it just like a raw pointer:
    ```cpp
    meshPtr->compute_voronoi_area();  // Use the object like a regular pointer
    ```
  - Alternatively, you can dereference using the `*` operator:
    ```cpp
    (*meshPtr).compute_voronoi_area();
    ```
    Basically the `*` gives you the object that the pointer points to.

- **Accessing Members**:
  - Smart pointers allow access to object members directly through the `->` operator, just like raw pointers:
    ```cpp
    HalfedgePtr nexthedge = current_hedge->getNextHalfEdge();
    ```

- **Moving Ownership**:
  - When you need to transfer ownership of an object managed by `std::unique_ptr`, you do so using `std::move()`:
    ```cpp
    std::unique_ptr<Mesh> mesh2 = std::move(meshPtr);  // meshPtr is now empty
    ```

- **Shared Access**:
  - With `std::shared_ptr`, you can have multiple pointers to the same object, and you can copy the pointer freely:
    ```cpp
    std::shared_ptr<Mesh> meshCopy = meshPtr;  // Now both meshPtr and meshCopy manage the same Mesh
    ```


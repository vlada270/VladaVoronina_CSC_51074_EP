#ifndef MESHPARTS_H
#define MESHPARTS_H
#include "Eigen/Core"
#include <memory>
#include <vector>
#include <iostream>

namespace MeshParts{


    //Forward declarations

    struct HalfEdge;
    struct PrimalFace;
    struct Vertex;

    typedef std::shared_ptr<PrimalFace> PrimalFacePtr;
	typedef std::shared_ptr<Vertex> VertexPtr;
	typedef std::shared_ptr<HalfEdge> HalfEdgePtr;


    struct Vertex {

        Vertex() = default;
        ~Vertex() = default;

        Vertex(int index): index(index){};

        int index;
        Eigen::Vector3d position;
        
        std::vector<std::weak_ptr<PrimalFace>> one_ring_faces;

        std::vector<std::weak_ptr<HalfEdge>> incoming_half_edges;
        double angle_defect = 0.;

        double gaussian_curvature = 0.;
        double voronoi_area = 0.;
        Eigen::Vector3d normal;

        std::vector<PrimalFacePtr> getOneRingFaces() const
		{
			std::vector<PrimalFacePtr> faces;

			for (auto face : one_ring_faces)
			{
				if (!face.expired())
				{
					faces.push_back(face.lock());
				}
				else
				{
					throw std::runtime_error("Face is deleted");
				}
			}
			return faces;
		}

        std::vector<HalfEdgePtr> getIncomingHalfEdges() const
		{
			std::vector<HalfEdgePtr> edges;
			for (auto edge : incoming_half_edges)
			{
				if (!edge.expired())
				{
					edges.push_back(edge.lock());
				}
				else
				{
					throw std::runtime_error("Edge is deleted");
				}
			}
			return edges;
		}

        void addOneRingFace(PrimalFacePtr one_ring_face)
		{
			this->one_ring_faces.push_back(one_ring_face);
		}

		void addIncomingHalfEdge(HalfEdgePtr incoming_half_edge)
		{
			this->incoming_half_edges.push_back(incoming_half_edge);
		}

		bool operator==(const Vertex &other) const;
		
        void print();

    };

    struct HalfEdge {

        HalfEdge() = default;
        ~HalfEdge() = default;

        HalfEdge(int index):index(index){};

        bool boundary = false;

        int index;
        // int sign_edge = 1;
        
        std::weak_ptr<Vertex> start;
		std::weak_ptr<Vertex> end;
		std::weak_ptr<HalfEdge> flip;
		std::weak_ptr<HalfEdge> next;
		std::weak_ptr<PrimalFace> primal_face;
		

        bool use_frame = false; //flag that indicates whether the frame field associated with the current face should be used for the integral.
        double hinge_connection_angle=0.;
        

        VertexPtr getStartVertex() const
		{
			if (!start.expired())
			{
				return start.lock();
			}
			else
			{
				throw std::runtime_error("Start Vertex is deleted");
			}
		}

		VertexPtr getEndVertex() const
		{
			if (!end.expired())
			{
				return end.lock();
			}
			else
			{
				throw std::runtime_error("End Vertex is deleted");
			}
		}

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

        void setStartVertex(VertexPtr start)
		{
			this->start = start;
		}

        void setEndVertex(VertexPtr end)
		{
			this->end = end;
		}

		void setFlipHalfEdge(HalfEdgePtr flip)
		{
			this->flip = flip;
		}

		void setNextHalfEdge(HalfEdgePtr next)
		{
			this->next = next;
		}

		void setPrimalFace(PrimalFacePtr primal_face)
		{
			this->primal_face = primal_face;
		}

        void print();

        bool operator==(const HalfEdge& other) const;
    };

    struct PrimalFace {
        PrimalFace () = default;
        // explicit DualFace(MatrixXd indices_face) : indices_face(indices_face) {}
        ~PrimalFace() = default;
        int index = -1;

        std::vector<std::weak_ptr<Vertex>> vertices_face;
		std::vector<std::weak_ptr<HalfEdge>> hedges_face;

        std::vector<int> indices_vertices;
        std::vector<int> indices_hedges;
        
        Eigen::Vector3d barycenter;
        Eigen::Vector3d circumcenter;
        
        Eigen::Vector3d normal;

        //ONLY CHANGE WRT TD4
        // frame field that is stored as a 2 x 3 matrix , i.e for every row we store an R^3 vector with a basis.
        Eigen::MatrixXd frame;
        
        std::vector<VertexPtr> getVertices() const
		{
			std::vector<VertexPtr> vertices;
			for (auto vertex : vertices_face)
			{
				if (!vertex.expired())
				{
					vertices.push_back(vertex.lock());
				}
				else
				{
					throw std::runtime_error("Vertex is deleted");
				}
			}
			return vertices;
		}

		std::vector<HalfEdgePtr> getHalfEdges() const
		{
			std::vector<HalfEdgePtr> hedges;
			for (auto hedge : hedges_face)
			{
				if (!hedge.expired())
				{
					hedges.push_back(hedge.lock());
				}
				else
				{
					throw std::runtime_error("Half Edge is deleted");
				}
			}
			return hedges;
		}

        void addVertex(VertexPtr vertex)
		{
			this->vertices_face.push_back(vertex);
		}

		void addHalfEdge(HalfEdgePtr hedge)
		{
			this->hedges_face.push_back(hedge);
		}

		bool operator==(const PrimalFace &other) const;

		void print();
    };

};

#endif
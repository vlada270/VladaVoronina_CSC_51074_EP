#include "MeshParts.h"

void MeshParts::Vertex::print() {
  std::cout << "position of vertex" << this->index << " is \n"<< this->position << std::endl;
}

void MeshParts::HalfEdge::print() {
  std::cout << "start: " << this->getStartVertex()->index<< " end: " << this->getEndVertex()->index << std::endl;
}

void MeshParts::PrimalFace::print()
{
    std::cout<<"vertices of face "<<this->index<<" are "<<this->indices_vertices[0]<<" "<<this->indices_vertices[1]<<" "<<this->indices_vertices[2]<<std::endl;
}

bool MeshParts::HalfEdge::operator==(const HalfEdge &other) const {
  return (
          (this->index == other.index) &&
          (this->getStartVertex() == other.getEndVertex()) &&
          (this->getFlipHalfEdge()->index == other.getFlipHalfEdge()->index) &&
          (this->getNextHalfEdge()->index == other.getNextHalfEdge()->index) &&
          (this->getPrimalFace()->index == other.getPrimalFace()->index) 
        );
}

bool MeshParts::Vertex::operator==(const Vertex &other) const {
  return (((position == other.position) && (index == other.index)));
}

bool MeshParts::PrimalFace::operator==(const PrimalFace &other) const {
  return ((indices_vertices == other.indices_vertices) &&
          (this->index == other.index) &&
          (this->normal == other.normal) &&
          (this->barycenter == other.barycenter) &&
          (this->indices_hedges == other.indices_hedges));
}
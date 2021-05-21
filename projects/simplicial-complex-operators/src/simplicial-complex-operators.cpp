// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        geometry->vertexIndices[i] = i;
    }

    // This is also a valid way to access indices (directly by mesh element).
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        geometry->vertexIndices[v] = idx;
        idx++;
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        geometry->edgeIndices[i] = i;
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        geometry->faceIndices[i] = i;
    }


    // You can more easily get the indices of mesh elements using the function getIndex(), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    std::vector<Eigen::Triplet<int>> triplets;

    for (Edge e : mesh->edges()) {
        for (Vertex v_i : e.adjacentVertices()) {
            triplets.push_back(Eigen::Triplet<int>(e.getIndex(), v_i.getIndex(), 1));
        }

    }
    SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());
    mat.setFromTriplets(triplets.begin(), triplets.end());

    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    std::vector<Eigen::Triplet<int>> triplets;

    for (Face f : mesh->faces()) {
        for (Edge e_i : f.adjacentEdges()) {
            triplets.push_back(Eigen::Triplet<int>(f.getIndex(), e_i.getIndex(), 1));
        }

    }
    SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    mat.setFromTriplets(triplets.begin(), triplets.end());

    return mat;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    Vector<size_t> verts = Vector<size_t>::Zero(mesh->nVertices(),1); // Column-Vector with nVertice rows
    for (size_t v : subset.vertices) {
        verts[v] = 1;
    }

    return verts;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> edges_v = Vector<size_t>::Zero(mesh->nEdges(), 1); // Column-Vector with nVertice rows
    for (size_t e : subset.edges) {
        edges_v[e] = 1;
    }

    return edges_v;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> face_v = Vector<size_t>::Zero(mesh->nFaces(), 1); // Column-Vector with nVertice rows
    for (size_t f : subset.faces) {
        face_v[f] = 1;
    }

    return face_v;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    MeshSubset result = subset.deepCopy();

    for (size_t v : subset.vertices) {
        // iterate over sparse matrix and add 1's
        for (SparseMatrix<size_t>::InnerIterator it(A0, v); it; ++it) {
            size_t e = it.row();
            result.addEdge(e);
        }
    }

    for (size_t e : result.edges) {
        for (SparseMatrix<size_t>::InnerIterator it(A1, e); it; ++it) {
            size_t f = it.row();
            result.addFace(f);
        }
    }

    // TODO
    return result; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}
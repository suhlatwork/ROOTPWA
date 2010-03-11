///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev::                             $: revision of last commit
// $Author::                          $: author of last commit
// $Date::                            $: date of last commit
//
// Description:
//      container class that holds all external information for
//      amplitude calculation
//      internally the decay process is represented as a graph using
//      the Boost Graph Library
//
//
// Author List:
//      Boris Grube          TUM            (original author)
//
//
//-------------------------------------------------------------------------


#include <list>

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/graphml.hpp>

#include "utilities.h"
#include "decayTopology.h"

	
using namespace std;
using namespace boost;
using namespace rpwa;


bool decayTopology::_debug = false;


decayTopology::decayTopology()
{ }


decayTopology::decayTopology(const decayTopology& topo)
{
  *this = topo;
}


decayTopology::decayTopology(const vector<particle*>&          fsParticles,
			     const vector<interactionVertex*>& vertices)
{
  constructDecay(fsParticles, vertices);
}


decayTopology::~decayTopology()
{ }


decayTopology&
decayTopology::operator = (const decayTopology& topo)
{
  if (this != &topo) {
    _graph           = topo._graph;
    _nodeProp        = topo._nodeProp;
    _edgeProp        = topo._edgeProp;
    _vertices        = topo._vertices;
    _fsParticles     = topo._fsParticles;
    _vertexNodeMap   = topo._vertexNodeMap;
    _particleEdgeMap = topo._particleEdgeMap;
  }
  return *this;
}


decayTopology&
decayTopology::constructDecay(const vector<particle*>&          fsParticles,
			      const vector<interactionVertex*>& vertices)
{
  clear();
  const unsigned int nmbVertices    = vertices.size();
  const unsigned int nmbFsParticles = fsParticles.size();
  if (nmbFsParticles == 0) {
    printWarn << "cannot construct decay topology without final state particles" << endl;
    clear();
    return *this;
  }
  if (nmbVertices == 0) {
    printWarn << "cannot construct decay topology without vertices" << endl;
    clear();
    return *this;
  }
  // create graph nodes for interaction vertices and store pointers
  _nodeProp = get(vertex_name, _graph);
  for (unsigned int i = 0; i < nmbVertices; ++i) {
    nodeDesc node = add_vertex(_graph);
    _nodeProp[node]             = vertices[i];
    _vertexNodeMap[vertices[i]] = node;
  }
  // create final state nodes
  for (unsigned int i = 0; i < nmbFsParticles; ++i) {
    nodeDesc node   = add_vertex(_graph);
    _nodeProp[node] = 0;
  }
  // create edges that connect the intercation nodes and store pointers to particles
  _edgeProp = get(edge_name, _graph);
  for (vertexNodeMapIt iFromVert = _vertexNodeMap.begin();
       iFromVert != _vertexNodeMap.end(); ++iFromVert)
    for (vertexNodeMapIt iToVert = _vertexNodeMap.begin();
	 iToVert != _vertexNodeMap.end(); ++iToVert) {
      interactionVertex* fromVertex = iFromVert->first;
      for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart) {
	interactionVertex* toVertex = iToVert->first;
	for (unsigned int iInPart = 0; iInPart < toVertex->nmbInParticles(); ++iInPart) {
	  particle* p = fromVertex->outParticles()[iOutPart];
	  if (p == toVertex->inParticles()[iInPart]) {
	    bool     inserted;
	    edgeDesc edge;
	    tie(edge, inserted) = add_edge(iFromVert->second, iToVert->second, _graph);
	    if (inserted) {
	      _edgeProp[edge]     = p;
	      _particleEdgeMap[p] = edge;
	    }
	  }
	}
      }
    }
  // create edges for final state particles and store pointers
  for (vertexNodeMapIt iFromVert = _vertexNodeMap.begin();
       iFromVert != _vertexNodeMap.end(); ++iFromVert) {
    interactionVertex* fromVertex = iFromVert->first;
    for (unsigned int iOutPart = 0; iOutPart < fromVertex->nmbOutParticles(); ++iOutPart)
      for (unsigned int iFsPart = 0; iFsPart < nmbFsParticles; ++iFsPart) {
	particle* p = fromVertex->outParticles()[iOutPart];
	if (p == fsParticles[iFsPart]) {
	  bool     inserted;
	  edgeDesc edge;
	  tie(edge, inserted) = add_edge(iFromVert->second, nmbVertices + iFsPart, _graph);
	  if (inserted) {
	    _edgeProp[edge]     = p;
	    _particleEdgeMap[p] = edge;
	  }
	}
      }
  }
  // sort nodes
  list<nodeDesc> sortedNodes;
  topological_sort(_graph, front_inserter(sortedNodes));
  for (list<nodeDesc>::iterator iNode = sortedNodes.begin();
       iNode != sortedNodes.end(); ++iNode) {
    interactionVertex* vertex = _nodeProp[*iNode];
    if (vertex)
      _vertices.push_back(vertex);
  }
  // copy final state particles
  _fsParticles = fsParticles;
  return *this;
}


bool
decayTopology::verifyTopology() const
{
  // check that decay topology is a tree of isobar decays
  bool          topologyOkay = true;
  nodeIndexType nodeIndex    = get(vertex_index, _graph);
  nodeIterator  iNode, iNodeEnd;
  unsigned int  countNodesWithNoInputEdge = 0;
  for (tie(iNode, iNodeEnd) = boost::vertices(_graph); iNode != iNodeEnd; ++iNode) {
    const unsigned int i = get(nodeIndex, *iNode);
    // check that each node has exactly 1 incoming edge (isobar)
    unsigned int nmbInEdges = in_degree(*iNode, _graph);
    if (nmbInEdges == 0) {
      ++countNodesWithNoInputEdge;
      if (_debug)
	printInfo << "assuming node[" << i << "] is production node" << endl;
    } else if (nmbInEdges != 1) {
      printWarn << "number of input edges of node[" << i << "] is "
    		<< nmbInEdges << " != 1" << endl;
      topologyOkay = false;
    } else if (_debug)
      printInfo << "number of input edges of node[" << i << "] is correct" << endl;
    if (countNodesWithNoInputEdge > 1) {
      printWarn << "there are " << countNodesWithNoInputEdge
		<< " nodes with no no input edges." << endl;
      topologyOkay = false;
    }
    // check that for each node the number of outgoing edges is either 2 (decay node) or 0 (final state node)
    unsigned int nmbOutEdges = out_degree(*iNode, _graph);
    if (nmbOutEdges == 0) {
      if (_nodeProp[*iNode]) {
	printWarn << "node[" << i << "] has no outgoing edges, "
		  << "but has a interactionVertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "final state node[" << i << "] is correct" << endl;
    } else if (nmbOutEdges == 2) {
      if (!_nodeProp[*iNode]) {
	printWarn << "node[" << i << "] has 2 outgoing edges, "
		  << "but no interactionVertex pointer assigned." << endl;
	topologyOkay = false;
      } else if (_debug)
	printInfo << "interaction node[" << i << "] is correct" << endl;
    } else {
      printWarn << "number of output edges of node[" << i << "] is "
		<< nmbOutEdges << " != 0 or 2" << endl;
      topologyOkay = false;
    }
  }
  return topologyOkay;
}


ostream&
decayTopology::print(ostream& out) const
{
  out << "decay topology nodes:" << endl;
  list<nodeDesc> sortedNodes;
  topological_sort(_graph, front_inserter(sortedNodes));
  for (list<nodeDesc>::const_iterator iNode = sortedNodes.begin();
       iNode != sortedNodes.end(); ++iNode) {
    const interactionVertex* vertex = _nodeProp[*iNode];
    if (vertex)
      out << "    decay ";
    else
      out << "    final state ";
    out << "node[" << *iNode << "] ";
    if (vertex)
      out << *vertex;
    else
      out << endl;
  }
  out << "decay topology particles:" << endl;
  edgeIterator iEdge, iEdgeEnd;
  for (tie(iEdge, iEdgeEnd) = edges(_graph); iEdge != iEdgeEnd; ++iEdge) {
    const particle* p = _edgeProp[*iEdge];
    out << "    particle" << *iEdge << " ";
    if (p) {
      out << p->name() << sign(p->charge()) << " ";
      if (!_nodeProp[target(*iEdge, _graph)])
	out << "(final state)";
      out << endl;
    } else
      out << "error! zero pointer to particle" << endl;
  }
  return out;
}


ostream&
decayTopology::writeGraphViz(ostream& out) const
{
  write_graphviz(out, _graph);
  return out;
}


void
decayTopology::clear()
{
  _graph.clear();
  _nodeProp = get(vertex_name, _graph);
  _edgeProp = get(edge_name, _graph);
  _vertices.clear();
  _fsParticles.clear();
  _vertexNodeMap.clear();
  _particleEdgeMap.clear();
}

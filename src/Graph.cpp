/*
 * teexGraph --- by Frank Takes --- https://github.com/franktakes/teexgraph
 * 
 * Graph class functions
 */

#include "Graph.h"
#include <Rcpp.h>
using namespace Rcpp;

// load the graph based on two vectors of equal length containing source and target
bool Graph::loadFromSourceTarget(std::vector<nodeidtype> & source, std::vector<nodeidtype> & target, const int until = -1) {
	int end = until;
	if(until == -1)
		end = (signed)source.size();
	for(int i=0; i<end; i++)
		addEdge(mapNode(source[i]), mapNode(target[i]));
    makeUndirected();		
	return true;
} // loadFromSourceTarget


// compute sixcycles (triangles in one-mode) and return them as a vector of paths
vector< vector<nodeidtype> > Graph::sixCycles(const nodeidtype from, const nodeidtype to) {
	vector< vector<nodeidtype> > results;
	int a = mapNode(from);
	int b = mapNode(to);
    int x,w,y,z;
    
    z = (signed)E[a].size();	
	for(int i=0; i<z; i++) { // neighbors of a, distance 1
		w = (signed)E[E[a][i]].size();
		for(int j=0; j<w; j++) { // distance 2 of a
			if(E[E[a][i]][j] != a) { 
				y = (signed)E[E[E[a][i]][j]].size();
				for(int k=0; k<y; k++) { // distance 3 of a
					if(E[E[E[a][i]][j]][k] != E[a][i]) { 
						x = (signed)E[E[E[E[a][i]][j]][k]].size();
						for(int l=0; l<x; l++) { // distance 4 of a
							if(E[E[E[E[a][i]][j]][k]][l] != a && E[E[E[E[a][i]][j]][k]][l] != E[E[a][i]][j] && edge(E[E[E[E[a][i]][j]][k]][l], b)) {		
								vector<nodeidtype> result;
								result.clear();
								result.push_back(from);
								result.push_back(revMapNode(E[a][i]));
								result.push_back(revMapNode(E[E[a][i]][j]));
								result.push_back(revMapNode(E[E[E[a][i]][j]][k]));	
								result.push_back(revMapNode(E[E[E[E[a][i]][j]][k]][l]));															
								result.push_back(to);	
								results.push_back(result);	
							} // if
						} // for l/x
					} // if	
				} // for k/y									
			} // if
		} // for j/s
	} // for i/z
	return results;
} // sixCycles

// compute fourcycles (squares) and return them as a vector of paths
vector< vector<nodeidtype> > Graph::fourCycles(const nodeidtype from, const nodeidtype to) {
	vector< vector<nodeidtype> > results;
	int a = mapNode(from);
	int b = mapNode(to);
    int w,z;
    
    z = (signed)E[a].size();	
	for(int i=0; i<z; i++) { // neighbors of a
		w = (signed)E[E[a][i]].size();
		for(int j=0; j<w; j++) { // distance 2 of a
			if(E[E[a][i]][j] != a && edge(E[E[a][i]][j], b)) {
				vector<nodeidtype> result;
				result.clear();
				result.push_back(from);
				result.push_back(revMapNode(E[a][i]));
				result.push_back(revMapNode(E[E[a][i]][j]));
				result.push_back(to);	
				results.push_back(result);											
			} // if
		} // for
	} // for
	return results;
} // fourCycles

// initialization of Graph object
Graph::Graph() {
    maxn = MAXN;
    clear();
} // Graph constructor

// initialization of Graph object
Graph::Graph(const int nmax) {
    assert(nmax > 0);
    maxn = nmax;
    clear();
} // Graph constructor

// erase the current Graph object
void Graph::clear() {
    nodeMapping.clear();
    revMapping.clear();
    E.assign(maxn, vector<int>(0));
    rE.assign(maxn, vector<int>(0));
    hasSelfLoop.assign(maxn, false);
    n = m = selfm = nexti = 0;
    loaded = sortedandunique = undirected = doneWCC = doneSCC = false;
} // clear


// map input file node node number to id in range [0,n-1]
int Graph::mapNode(const nodeidtype i) {
    if(nodeMapping.find(i) == nodeMapping.end()) {
        nodeMapping[i] = nexti;
        revMapping[nexti] = i;
        nexti++;
    } // if
    return nodeMapping[i];
} // mapNode


// reverse map node id to original node number
nodeidtype Graph::revMapNode(const int i) {
    return revMapping[i];
} // revMapNode


// add an edge to the graph and update the necessary statistics
bool Graph::addEdge(const int u, const int v) {

    // check if node is within bounds set in Graph.h
    if(u < 0 || u >= maxn || v < 0 || v >= maxn) {
        return false;
    }

    // note if it is a self-loop
    if(u == v) {
        selfm++;
        hasSelfLoop[u] = true;
    }

    // increment edge count
    m++;

    // increment node count if one of the two nodes is first seen
    if(rE[v].size() == 0 && E[v].size() == 0) {
        n++;
    }
    if(v != u && rE[u].size() == 0 && E[u].size() == 0) {
        n++;
    }

    E[u].push_back(v);
    rE[v].push_back(u);

    sortedandunique = false;
    undirected = false;
    doneWCC = false;
    doneSCC = false;
    return true;
} // addEdge


// load an undirected graph from a file in edge list format: [u v]
bool Graph::loadUndirected(const string filename) {
    if(!loadDirected(filename))
        return false;
    makeUndirected();
    return true;
} // loadUndirected

// make graph undirected (i.e., if (u,v) and !(v,u), then add (u,v))
void Graph::makeUndirected() {
    int z, oldm = m;

    if(undirected) {
        clog << "Graph is already undirected." << endl;
        return;
    }

    //clog << "Making graph undirected (m = " << m << ")..." << endl;

    // add all links
    sortedandunique = false;

    for(int i = 0; i < n; i++) {
        z = E[i].size();
        for(int j = 0; j < z; j++)
            E[E[i][j]].push_back(i);
        rE[i].clear();
    } // for

    sortEdgeList(); // needed to remove duplicates introduced in previous step
    doneWCC = doneSCC = false;
    undirected = true;
    if(m != oldm && m != oldm * 2) {
        //cerr << "  WARNING: number of edges is not equal to (twice the) number of input lines."
        //        << endl << "  Verify that the graph is actually undirected." << endl;
    }
    //clog << "Undirected-making done (m = " << m << ")." << endl;
} // makeUndirected

// load a graph from a file in edge list format: [u v]
bool Graph::loadDirected(const string filename) {
    nodeidtype u, v;
    long edgesAdded = 0, edgesSkipped = 0;
    ifstream fin;

    clog << endl << "Loading graph from " << filename << " ..." << endl;

    // check if not already loaded
    if(loaded) {
        cerr << "Error: a graph is already loaded. Clear it first." << endl;
        return false;
    }

    // load the file
    fin.open(filename.c_str());
    if(!fin.is_open()) {
        cerr << "Error: file " << filename << " not found." << endl;
        return false;
    }

    // ignore first lines that do not start with an integer character 
    char c = fin.peek();
    while(!(
            (c >= '0' && c <= '9') ||
            (c >= 'a' && c <= 'z') ||
            (c >= 'A' && c <= 'Z')
            )) {
        while(c != '\n')
            c = fin.get();
        c = fin.peek();
    }

    // load the edge list
    while(fin >> u >> v) {
        if(m % 10000000 == 0 && m != 0)
            clog << "   - " << m << " edges loaded so far..." << endl;
        if(addEdge(mapNode(u), mapNode(v)))
            edgesAdded++;
        else
            edgesSkipped++;
    }

    clog << "- " << edgesAdded << " edges added (m = " << m << ") in total\n- "
            << edgesSkipped << " edges skipped" << endl;
    if(edgesSkipped - selfm > 0)
        clog << " (out-of-bounds, increase maxn in Graph.h!)";
    clog << "- " << selfm << " self-edges added" << endl;
    clog << endl;

    loaded = true;

    // succesful if we didnt have to skip edges	
    if(edgesSkipped == 0) {
        sortEdgeList();
        clog << "Loading done." << endl << endl;
        return true;
    }
    cerr << "Loading failed." << endl << endl;
    clear();
    return false;
} // loadDirected


// check if there is an edge from a to b - O(log(outdegree(a)))
bool Graph::edge(const int a, const int b) {
    if(!sortedandunique)
        sortEdgeList();
    int first = 0, last = E[a].size() - 1, mid;
    while(first <= last) {
        mid = (first + last) / 2;
        if(b > E[a][mid])
            first = mid + 1;
        else if(b < E[a][mid])
            last = mid - 1;
        else
            return true;
    }
    return false;
} // edge


// check if there is an edge from a to b without requiring sortid list - O(outdegree(a))
bool Graph::edgeSlow(const int a, const int b) {
    const int z = E[a].size();
    for(int j = 0; j < z; j++)
        if(E[a][j] == b)
            return true;
    return false;
} // edgeSlow


// sort edge list so that O(log(outdegree(a))) queries edge(a,b) are possible
void Graph::sortEdgeList() {
    //clog << "Sorting edge list..." << endl;
    int removed = 0, begin;
    if(!sortedandunique) {
        m = 0;
        // TODO: make parallel, atomic m and removed update
        for(int i = 0; i < n; i++) {
            sort(E[i].begin(), E[i].end());
            begin = E[i].size();
            E[i].erase(unique(E[i].begin(), E[i].end()), E[i].end());
            sort(rE[i].begin(), rE[i].end());
            rE[i].erase(unique(rE[i].begin(), rE[i].end()), rE[i].end());
            m += (long) E[i].size();
            removed += (E[i].size() - begin);
        }

        sortedandunique = true;
        if(removed > 0)
            clog << "- Total of " << removed << " parallel edges have been removed." << endl;
        //clog << "Sorting done." << endl;
    } else
        clog << "Was already done." << endl;
} // sortEdgeList



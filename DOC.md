*To tommorow:*     
    I ve written the rabbani method's call at the `.cpp` file where the `main()` function exists.  
    We should keep track of the methods' calling, and implement those method ONE BY ONE.  
    DONT care about others. COOL.  
    Chill bro, ... Chill, ...  
    Success Never Ever comes EASY.  
    Road to Top Tier AI conference...  


# Flow of H-Region Merging 

- `OpenLAS()`: \{  
	create `HRM3D` passing the data `path`  
	`run()` method `HRM3D` with params $[scale\_param, alpha, search\_radius, region\_size, mode]$  
 \}

- `HRM3D` instantiation:  
	two functions. cloud. 

- `HRM3D` method `run()`:  
	create `RAG_Cloud` with params $[cloud, alpha, search\_radius, region\_size, mode]$  
    call `RAG_Cloud`'s `run()` method with params $[scale\_params, mode, region\_size]$

`RAG_Cloud`'s `run()` method. Params: $[thresh, mode, region\_size]$



# TODO

## Main stream

1. [x] Quick: copy pasting the `HRM3D`'s Constructor for `RabbaniRM3D`. 
2. [x] Quick: implement the `run()` method for Rab`baniRM3D`. 
3. [x] Hard:  Copy pasting `RAG_Cloud` Constructor for `RabbaniRAG`. 
4. [x] Hard:  Tailor the class for *Rabbani Criteria*.
5. [ ] Hard:  Debug and Run. 

## Added during the Process

1. [x] Quick: Implementing the `connect_vertex_and_edge()` method in the `RabbaniRAG` class. 
2. [x] Hard:  Implementing the `initial_vertice()` inside `RabbaniRAG` class working as a helper. 
   1. [ ] <span style="color:red">CHECK: `v->normal = n;`, I am not sure whether this is op on `RabbaniVertex`.</span> Should come back later ...
3. [x] Hard: Implementing the `kd_tree_radius_consturct_topo(search_radius)` method in the `RabbaniRAG` class.
4. [x] Hard: Implementing the `init_edges(edge_count)` method in `RabbaniRAG` class.  
   1. [x] <span style="color:red">Discover the function of `merging_criteria(f, f, edge)`, returning the new weight.
    Which in this case, should be the angle difference. </span>
5. [x] Hard: Implementing the edge class `RabbaniEdge`. 
6. [x] Hard: Implementing the function helper - `updateFeatures()` & `mergingCriteria()`

# DOC

1. Implemented the `OpenLASRabbani()` with generic flow. 
2. Implemneted the `RabbaniRAG` class with basic Constructor which is used in `RabbaniRM3D`'s `run()` method. 
3. Adding a `RAGCallback.h` file for efficient use of helper functions.
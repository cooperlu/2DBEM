% function demo0
%DEMO0 a very simple example to start with -- mesh a square
%domain with a square hold cut from its centre.

    fprintf(1, [ ...
' A very simple example to start with -- construct a mesh for \n', ...
' a simple square domain with a square hole cut from its cen- \n', ...
' tre. The geometry is specified as a Planar Straight-Line \n', ...
' Graph (PSLG) -- a list of xy coordinates, or "nodes", and a \n', ...
' list of straight-line connections between nodes, or "edges".\n', ...
' The REFINE2 routine is used to build a triangulation of the \n', ...
' domain that: (a) conforms to the geometry, and (b) contains \n', ...
' only "nicely" shaped triangles. In the second panel, a mesh \n', ...
' that additionally satisfies "mesh-size" constrains is cons- \n', ...
' structed -- '
        ] ) ;

%------------------------------------------- setup geometry
    
    node = [                % list of xy "node" coordinates
        0, 0                % outer square
        9, 0
        9, 9
        0, 9 
        4, 4                % inner square
        5, 4
        5, 5
        4, 5 ] ;
    
    edge = [                % list of "edges" between nodes
        1, 2                % outer square 
        2, 3
        3, 4
        4, 1 
        5, 6                % inner square
        6, 7
        7, 8
        8, 5 ] ;

%------------------------------------------- call mesh-gen.
   [vert,etri, ...
    tria,tnum] = refine2(node,edge) ;

%------------------------------------------- draw tria-mesh
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    
%------------------------------------------- call mesh-gen.
    hfun = +.5 ;            % uniform "target" edge-lengths

   [vert,etri, ...
    tria,tnum] = refine2(node,edge,[],[],hfun) ;

%------------------------------------------- draw tria-mesh
    figure;
    patch('faces',tria(:,1:3),'vertices',vert, ...
        'facecolor','w', ...
        'edgecolor',[.2,.2,.2]) ;
    hold on; axis image off;
    patch('faces',edge(:,1:2),'vertices',node, ...
        'facecolor','w', ...
        'edgecolor',[.1,.1,.1], ...
        'linewidth',1.5) ;
    
    drawnow;
    
    set(figure(1),'units','normalized', ...
        'position',[.05,.50,.30,.35]) ;
    set(figure(2),'units','normalized', ...
        'position',[.35,.50,.30,.35]) ;
    
% end
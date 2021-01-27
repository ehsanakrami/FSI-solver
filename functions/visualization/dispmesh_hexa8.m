function dispmesh_hexa8(T_E,T_X,colorstr,alpha)

Ne = size(T_E,1) ;

Xface = zeros(4,6*Ne) ;
Yface = zeros(4,6*Ne) ;
Zface = zeros(4,6*Ne) ;

for ee = 1:Ne
    
    nodes     = T_E(ee,:)   ;
    
    Xface_elt = [T_X(nodes([1 2 3 4]),1) , T_X(nodes([1 4 8 5]),1) , T_X(nodes([1 2 6 5]),1) , T_X(nodes([2 3 7 6]),1) , T_X(nodes([3 4 8 7]),1) , T_X(nodes([5 6 7 8]),1)] ;
    Yface_elt = [T_X(nodes([1 2 3 4]),2) , T_X(nodes([1 4 8 5]),2) , T_X(nodes([1 2 6 5]),2) , T_X(nodes([2 3 7 6]),2) , T_X(nodes([3 4 8 7]),2) , T_X(nodes([5 6 7 8]),2)] ;
    Zface_elt = [T_X(nodes([1 2 3 4]),3) , T_X(nodes([1 4 8 5]),3) , T_X(nodes([1 2 6 5]),3) , T_X(nodes([2 3 7 6]),3) , T_X(nodes([3 4 8 7]),3) , T_X(nodes([5 6 7 8]),3)] ;
    
    Xface(:,(6*ee-5):(6*ee)) = Xface_elt ;
    Yface(:,(6*ee-5):(6*ee)) = Yface_elt ;
    Zface(:,(6*ee-5):(6*ee)) = Zface_elt ;
    
end

patch(Xface,Yface,Zface,colorstr,'FaceAlpha',alpha) ;

axis equal
view(3)
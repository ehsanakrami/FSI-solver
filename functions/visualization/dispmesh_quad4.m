function dispmesh_quad4(T_E,T_X,colorstr,alpha)

% Affichage du maillage initial

Ne = size(T_E,1) ;

Xface = zeros(4,Ne) ;
Yface = zeros(4,Ne) ;
Zface = zeros(4,Ne) ;

for ee = 1:Ne

    nodes = T_E(ee,:)   ;

    Xface_elt = T_X(nodes([1 2 3 4]),1) ;
    Yface_elt = T_X(nodes([1 2 3 4]),2) ;
    Zface_elt = T_X(nodes([1 2 3 4]),3) ;

    Xface(:,ee) = Xface_elt ;
    Yface(:,ee) = Yface_elt ;
    Zface(:,ee) = Zface_elt ;

end

patch(Xface,Yface,Zface,colorstr,'FaceAlpha',alpha) ;

axis equal
view(3)
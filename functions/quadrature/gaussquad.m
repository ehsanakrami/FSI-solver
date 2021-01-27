function [wg,pg] = gaussquad(nn)

if nn == 8
    
    % 8 point integration
    
    a = 1/sqrt(3) ;
    
    wg = [ 1 ;
        1 ;
        1 ;
        1 ;
        1 ;
        1 ;
        1 ;
        1  ] ;
    
    pg = [ a , a , a    ;
        a , a , -a   ;
        a , -a , a   ;
        a , -a , -a  ;
        -a , a , a   ;
        -a , a , -a  ;
        -a , -a , a  ;
        -a , -a , -a   ] ;
    
end

if nn == 1
    
    % 1 point integration
    
    wg = 8 ;
    
    pg = [0 , 0 , 0] ;
    
end

if nn == 27
    
    % 27 point integration
    
    c = sqrt(3/5) ;
    
    w0 = 8/9 ;
    w  = 5/9 ;
    
    wg = [ w*w*w    ;
        w*w*w    ;
        w*w*w0   ;
        w*w*w    ;
        w*w*w    ;
        w*w*w0   ;
        w*w0*w   ;
        w*w0*w   ;
        w*w0*w0  ;
        w*w*w    ;
        w*w*w    ;
        w*w*w0   ;
        w*w*w    ;
        w*w*w    ;
        w*w*w0   ;
        w*w0*w   ;
        w*w0*w   ;
        w*w0*w0  ;
        w0*w*w   ;
        w0*w*w   ;
        w0*w*w0  ;
        w0*w*w   ;
        w0*w*w   ;
        w0*w*w0  ;
        w0*w0*w  ;
        w0*w0*w  ;
        w0*w0*w0  ] ;
    
    pg = [  c  ,  c  ,  c  ;
        c  ,  c  , -c  ;
        c  ,  c  ,  0  ;
        c  , -c  ,  c  ;
        c  , -c  , -c  ;
        c  , -c  ,  0  ;
        c  ,  0  ,  c  ;
        c  ,  0  , -c  ;
        c  ,  0  ,  0  ;
        -c  ,  c  ,  c  ;
        -c  ,  c  , -c  ;
        -c  ,  c  ,  0  ;
        -c  , -c  ,  c  ;
        -c  , -c  , -c  ;
        -c  , -c  ,  0  ;
        -c  ,  0  ,  c  ;
        -c  ,  0  , -c  ;
        -c  ,  0  ,  0  ;
        0  ,  c  ,  c  ;
        0  ,  c  , -c  ;
        0  ,  c  ,  0  ;
        0  , -c  ,  c  ;
        0  , -c  , -c  ;
        0  , -c  ,  0  ;
        0  ,  0  ,  c  ;
        0  ,  0  , -c  ;
        0  ,  0  ,  0   ] ;
    
end
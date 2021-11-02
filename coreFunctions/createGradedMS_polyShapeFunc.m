function [psFuncPos,psFuncNeg] = createGradedMS_polyShapeFunc(lType)
%WRITELATTICETOPOLOGY
if nargin <1
    lType = 8;
end
%To suppress warning when some voids vanish at some thickness
warning('off','all');
%nodes = [0 0;1 0;1 1;0 1;0.5 0;1 0.5;0.5 1;0 0.5;0.5 0.5];
switch lType
    case 1  %Box
        xs = @(x) [0;1;1;0];
        ys = @(x) [0;0;1;1];
        xs1 = @(x) [x;1-x;1-x;x];ys1 = @(x) [1-x;1-x;x;x];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape(xs1(x),ys1(x));
    case 2  %Plus
        xs = @(x) [0.5-x/2;0.5+x/2;0.5+x/2;1;1;0.5+x/2;0.5+x/2;0.5-x/2;0.5-x/2;0;0;0.5-x/2];
        ys = @(x) [0;0;0.5-x/2;0.5-x/2;0.5+x/2;0.5+x/2;1;1;0.5+x/2;0.5+x/2;0.5-x/2;0.5-x/2];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape();
    case 3  %inverted-Z
        a = @(x) sqrt(2)*x;
        xs = @(x) 1-[0;1;1   ;1-(a(x)-x);x+a(x);1     ;1;0;0     ;a(x)-x;1-(a(x)+x);0];
        ys = @(x) [0;0;a(x);x         ;x     ;1-a(x);1;1;1-a(x);1-x   ;1-x       ;a(x)];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape();
    case 4  %Z
        a = @(x) sqrt(2)*x;
        xs = @(x) [0;1;1   ;1-(a(x)-x);x+a(x);1     ;1;0;0     ;a(x)-x;1-(a(x)+x);0];
        ys = @(x) [0;0;a(x);x         ;x     ;1-a(x);1;1;1-a(x);1-x   ;1-x       ;a(x)];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape();
    case 5  %N
        a = @(x) sqrt(2)*x;
        ys = @(x) 1-[0;1;1   ;1-(a(x)-x);x+a(x);1     ;1;0;0     ;a(x)-x;1-(a(x)+x);0];
        xs = @(x) [0;0;a(x);x         ;x     ;1-a(x);1;1;1-a(x);1-x   ;1-x       ;a(x)];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape();
    case 6  %inverted-N
        a = @(x) sqrt(2)*x;
        ys = @(x) [0;1;1   ;1-(a(x)-x);x+a(x);1     ;1;0;0     ;a(x)-x;1-(a(x)+x);0];
        xs = @(x) [0;0;a(x);x         ;x     ;1-a(x);1;1;1-a(x);1-x   ;1-x       ;a(x)];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape();
    case 7  %X
        a = @(x) sqrt(2)*x;
        xs = @(x) [0;a(x);0.5;1-a(x);1;1;0.5+a(x);1;1;1-a(x);0.5;a(x);0;0;0.5-a(x);0];
        ys = @(x) [0;0;0.5-a(x);0;0;a(x);0.5;1-a(x);1;1;0.5+a(x);1;1;1-a(x);0.5;a(x)];
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape();
    case 8  %Hexagonal
        s = 1;lx = 3*s;ly = sqrt(3)*s;a = @(x) x/sqrt(3);
        xs = @(x) [0     ;s/2-a(x);s-2*a(x);2*s+2*a(x);5*s/2+a(x);lx    ;...
            lx    ;5*s/2+a(x);2*s+2*a(x);s-2*a(x);s/2-a(x);0     ]/ly;
        ys = @(x) [ly/2-x;ly/2-x  ;0       ;0         ;ly/2-x    ;ly/2-x;...
            ly/2+x;ly/2+x    ; ly      ;ly       ;ly/2+x  ;ly/2+x]/ly;
        xs1 = @(x) [s/2+2*a(x);s+a(x);2*s-a(x);5*s/2-2*a(x);2*s-a(x);s+a(x)]/ly;
        ys1 = @(x) [ly/2      ;  x   ; x      ;ly/2        ;ly-x    ;ly-x]/ly;
        %Division by ly to scale smallest dimension to unity
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape(xs1(x),ys1(x));
    case 9 %Kagome = isotropic
        s = 1;lx = 3*s;ly = sqrt(3)*s;xf1 = lx/4;
        a = @(x) x/sqrt(3);r = @(x) ly/lx*x;r1 = @(x) lx/ly*x;
        xs = @(x)[0     ;max(0,lx/4-x)  ;max(0,lx/4-x);min(lx,lx/4+x);min(lx,lx/4+x) ;max(0,lx/2-r1(x));...
            min(lx,lx/2+r1(x));max(0,3*lx/4-x);max(0,3*lx/4-x);min(lx,3*lx/4+x);min(lx,3*lx/4+x);lx    ;...
            lx     ;min(lx,3*lx/4+x)  ;min(lx,3*lx/4+x);max(0,3*lx/4-x);max(0,3*lx/4-x) ;min(lx,lx/2+r1(x));...
            max(0,lx/2-r1(x));min(lx,lx/4+x);lx/4+x;max(0,lx/4-x);max(0,lx/4-x);0]/ly;
        %right cavity
        ys = @(x)[ly/2-x;max(0,ly/4-x+r(x));0     ;0     ;max(0,ly/4-x-r(x));    0      ;...
            0      ;max(0,ly/4-x-r(x));0       ;0       ;max(0,ly/4-x+r(x));ly/2-x;...
            ly/2+x;min(ly,3*ly/4+x-r(x));ly    ;ly     ;min(ly,3*ly/4+x+r(x));    ly ;...
            ly      ;min(ly,3*ly/4+x+r(x));ly       ;ly       ;min(ly,3*ly/4+x-r(x));ly/2+x;]/ly;%right cavity
        xs1 = {    @(x) [min(lx,lx/4+x);lx/2;max(0,3*lx/4-x);max(0,3*lx/4-x);lx/2;min(lx,lx/4+x)]/ly;%Central cavity
            @(x) (2*x<lx/4-x)*[2*x;lx/4-x;lx/4-x]/ly;%left cavity
            @(x) (2*x<lx/4-x)*[lx-2*x;3*lx/4+x;3*lx/4+x]/ly};
        ys1 = { @(x)[ly/4+x-r(x);2*a(x);ly/4+x-r(x);3*ly/4-x+r(x);ly-2*a(x);3*ly/4-x+r(x)]/ly;%Central cavity
            @(x) (2*x<lx/4-x)*[ly/2;ly/4+x+r(x);3*ly/4-x-r(x)]/ly;%left cavity
            @(x) (2*x<lx/4-x)*[ly/2;3*ly/4-x-r(x);ly/4+x+r(x)]/ly};
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = cell(size(xs1,1),1);
        for i = 1:size(xs1,1)
            psFuncNeg{i} = @(x) polyshape(xs1{i}(x),ys1{i}(x));
        end
    case 10 %Diamond = orthotropic
        s = 1;lx = 3*s;ly = sqrt(3)*s;
        r = @(x) ly/lx*x;r1 = @(x) lx/ly*x;
        xs = @(x)[0;lx/2-r1(x);lx/2+r1(x);lx;lx ;lx/2+r1(x);lx/2-r1(x);0]/ly;
        
        ys = @(x)[ly/2-x;0;0;ly/2-x;ly/2+x;ly;ly;ly/2+x]/ly;
        xs1 = {    @(x) [2*x;lx/2-x;lx/2-x]/ly;%left cavity
            @(x) [lx-2*x;lx/2+x;lx/2+x]/ly};%right cavity
        ys1 = {@(x) [ly/2;x+r(x);ly-x-r(x)]/ly;%left cavity
            @(x) [ly/2;ly-x-r(x);x+r(x)]/ly};%right cavity
        %Division by ly to scale smallest dimension to unity
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = cell(size(xs1,1),1);
        for i = 1:size(xs1,1)
            psFuncNeg{i} = @(x) polyshape(xs1{i}(x),ys1{i}(x));
        end
    case 11 %Triangular honeycomb = almost isotropic
        s = 1;lx = 3*s;ly = sqrt(3)*s;
        r = @(x) ly/lx*x;r1 = @(x) lx/ly*x;a = @(x) x/sqrt(3);
        xs = @(x)[0;x;x;lx/2-r1(x);...
            lx/2+r1(x);lx-x;lx-x ;lx;...
            lx;lx-x;lx-x ;lx/2+r1(x);...
            lx/2-r1(x);x;x;0]/ly;
        
        ys = @(x)[0;0;ly/2-2*a(x)-r(x);0;...
            0;ly/2-2*a(x)-r(x);0;0;...
            ly;ly;ly/2+2*a(x)+r(x);ly;...
            ly;ly/2+2*a(x)+r(x);ly;ly]/ly;
        xs1 = {    @(x) [2*x;lx/2-x;lx/2-x]/ly;%left cavity
            @(x) [lx-2*x;lx/2+x;lx/2+x]/ly};%right cavity
        ys1 = {@(x) [ly/2;x+r(x);ly-x-r(x)]/ly;%left cavity
            @(x) [ly/2;ly-x-r(x);x+r(x)]/ly};%right cavity
        %Division by ly to scale smallest dimension to unity
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = cell(size(xs1,1),1);
        for i = 1:size(xs1,1)
            psFuncNeg{i} = @(x) polyshape(xs1{i}(x),ys1{i}(x));
        end
    case 12 %Diamond-Hexagonal = orthotropic
        s = 1;lx = 2*s;ly = sqrt(3)*s;a = @(x) x/sqrt(3);
        xs = @(x) [0     ;lx/2-s/2-a(x);lx/2-2*a(x);lx/2+2*a(x);lx/2+s/2+a(x);lx;...
            lx;lx/2+s/2+a(x);lx/2+2*a(x);lx/2-2*a(x);lx/2-s/2-a(x);0]/ly;
        ys = @(x) [ly/2-x;ly/2-x       ;0   ;0     ;ly/2-x    ;ly/2-x;...
            ly/2+x;ly/2+x   ; ly      ;ly       ;ly/2+x  ;ly/2+x]/ly;
        xs1 = @(x) [lx/2-s/2+2*a(x);lx/2;lx/2+s/2-2*a(x);lx/2]/ly;
        ys1 = @(x) [ly/2 ; 2*x ; ly/2 ;ly-2*x]/ly;
        %Division by ly to scale smallest dimension to unity
        psFuncPos = @(x) polyshape(xs(x),ys(x));
        psFuncNeg = @(x) polyshape(xs1(x),ys1(x));
    case 13 %Diamond-Octagon = isotropic orthotropic
        s = 1/(sqrt(2) + 1);x = (sqrt(2) + 1)*s;ly = (sqrt(2) + 1)*s;f = s/sqrt(2);
        %To be implemented
        psFuncPos = @(x) polyshape();
        psFuncNeg = @(x) polyshape();
    otherwise%To write a temporary lattice connectivity
        
end
end


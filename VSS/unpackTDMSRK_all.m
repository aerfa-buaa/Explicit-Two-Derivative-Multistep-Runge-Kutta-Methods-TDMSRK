function [A,Ahat,v,vhat,d,b,L] =  unpackTDMSRK_all(x, chars, w)
%% TDMSRK2249
L=1;
if chars==2249
    xx=[0.0776580383519972,0.301371373051056,0.199267993301564, 0.0907759259244949,0.0859182382169720];
    xx=   [0.0776580383519972,0.301371373051056,0.199267993301564,0.0907759259244949,0.0859182382169720,-0.946487048652487];
    b1 = xx(1);   v3 = xx(2); vv3 = xx(3); a21 = xx(4);  d31 = xx(5);  w1=w(1);
    % b1 = x(1);   v3 = x(2); vv3 = x(3); a21 = x(4);  a22 = x(5); d31 = x(6);  w1=w(1);
    
%     a21=0.0926*w1^(-2.352)    ;
%     b1=  0.0783*w1^(  -1.177)    ;
%     d31= 0.0586*w1^(-1.886  )+0.0273  ;
%     v3=-0.00174 *w1^(4.867  )+ 0.302  ;
%     vv3 = 0.0029*w1^(-3.627 )+ 0.196 ;
    
    a21= 0.07982*w1^(-2.688 ) + 0.01096 ;
    b1=  0.06514*w1^( -1.418  ) + 0.0125 ;
    d31= 0.06017*w1^( -1.835   ) + 0.02575  ;
    v3= 0.3985*exp(- 0.2794*w1) *w1^0.2483  ;
    vv3 =0.1651*exp( 0.1884*w1) *w1^(-0.2552)    ;
    
    a22= (-a21)*w1 + d31*w1 + (3*a21*w1^3 - d31*w1^3)^(1/3);
    
    b2=1 - b1 ;      d32=1 - d31;
    %     [1,0.0776580383519972,0.301371373051056,0.199267993301564,
    %      0.0907759259244949,0.566387706162208,0.0859182382169720,-0.946487048652487]
    %     [0.0776580383519972,0.301371373051056,0.199267993301564,
    %         0.0907759259244949,0.0859182382169720,-0.946487048652487]
    v1=(1/2).*w1.^(-4).*(1+(-4).*a22.^3.*v3+(-12).*a22.^2.*vv3+2.*w1+(-6) ...
        .*a22.^2.*v3.*w1+(-12).*a21.*a22.^2.*v3.*w1+12.*a22.^2.*d31.*v3.* ...
        w1+(-12).*a22.*vv3.*w1+(-24).*a21.*a22.*vv3.*w1+24.*a22.*d31.* ...
        vv3.*w1+(-12).*a21.*a22.*v3.*w1.^2+(-12).*a21.^2.*a22.*v3.*w1.^2+ ...
        12.*a22.*d31.*v3.*w1.^2+24.*a21.*a22.*d31.*v3.*w1.^2+(-12).*a22.* ...
        d31.^2.*v3.*w1.^2+(-12).*a21.*vv3.*w1.^2+(-12).*a21.^2.*vv3.* ...
        w1.^2+12.*d31.*vv3.*w1.^2+24.*a21.*d31.*vv3.*w1.^2+(-12).*d31.^2.* ...
        vv3.*w1.^2+(-6).*a21.^2.*v3.*w1.^3+(-4).*a21.^3.*v3.*w1.^3+12.* ...
        a21.*d31.*v3.*w1.^3+12.*a21.^2.*d31.*v3.*w1.^3+(-6).*d31.^2.*v3.* ...
        w1.^3+(-12).*a21.*d31.^2.*v3.*w1.^3+4.*d31.^3.*v3.*w1.^3+b1.*  w1.^4);
    v2=(1/2).*w1.^(-3).*((-1)+4.*a22.^3.*v3+12.*a22.^2.*vv3+(-2).*w1+6.* ...
        a22.^2.*v3.*w1+12.*a21.*a22.^2.*v3.*w1+(-12).*a22.^2.*d31.*v3.*w1+ ...
        12.*a22.*vv3.*w1+24.*a21.*a22.*vv3.*w1+(-24).*a22.*d31.*vv3.*w1+ ...
        12.*a21.*a22.*v3.*w1.^2+12.*a21.^2.*a22.*v3.*w1.^2+(-12).*a22.* ...
        d31.*v3.*w1.^2+(-24).*a21.*a22.*d31.*v3.*w1.^2+12.*a22.*d31.^2.* ...
        v3.*w1.^2+12.*a21.*vv3.*w1.^2+12.*a21.^2.*vv3.*w1.^2+(-12).*d31.* ...
        vv3.*w1.^2+(-24).*a21.*d31.*vv3.*w1.^2+12.*d31.^2.*vv3.*w1.^2+2.* ...
        w1.^3+(-2).*v3.*w1.^3+6.*a21.^2.*v3.*w1.^3+4.*a21.^3.*v3.*w1.^3+( ...
        -12).*a21.*d31.*v3.*w1.^3+(-12).*a21.^2.*d31.*v3.*w1.^3+6.* ...
        d31.^2.*v3.*w1.^3+12.*a21.*d31.^2.*v3.*w1.^3+(-4).*d31.^3.*v3.* w1.^3+b1.*w1.^4);
    vv1=(1/12).*w1.^(-4).*(3+(-12).*a22.^3.*v3+(-36).*a22.^2.*vv3+4.*w1+( ...
        -12).*a22.^2.*v3.*w1+(-36).*a21.*a22.^2.*v3.*w1+36.*a22.^2.*d31.* ...
        v3.*w1+(-24).*a22.*vv3.*w1+(-72).*a21.*a22.*vv3.*w1+72.*a22.*d31.* ...
        vv3.*w1+(-24).*a21.*a22.*v3.*w1.^2+(-36).*a21.^2.*a22.*v3.*w1.^2+ ...
        24.*a22.*d31.*v3.*w1.^2+72.*a21.*a22.*d31.*v3.*w1.^2+(-36).*a22.* ...
        d31.^2.*v3.*w1.^2+(-24).*a21.*vv3.*w1.^2+(-36).*a21.^2.*vv3.* ...
        w1.^2+24.*d31.*vv3.*w1.^2+72.*a21.*d31.*vv3.*w1.^2+(-36).*d31.^2.* ...
        vv3.*w1.^2+(-12).*a21.^2.*v3.*w1.^3+(-12).*a21.^3.*v3.*w1.^3+24.* ...
        a21.*d31.*v3.*w1.^3+36.*a21.^2.*d31.*v3.*w1.^3+(-12).*d31.^2.*v3.* ...
        w1.^3+(-36).*a21.*d31.^2.*v3.*w1.^3+12.*d31.^3.*v3.*w1.^3+b1.* w1.^4);
    vv2=(1/12).*w1.^(-2).*(3+(-12).*a22.^3.*v3+(-36).*a22.^2.*vv3+8.*w1+( ...
        -24).*a22.^2.*v3.*w1+(-36).*a21.*a22.^2.*v3.*w1+36.*a22.^2.*d31.* ...
        v3.*w1+(-48).*a22.*vv3.*w1+(-72).*a21.*a22.*vv3.*w1+72.*a22.*d31.* ...
        vv3.*w1+6.*w1.^2+(-12).*a22.*v3.*w1.^2+(-48).*a21.*a22.*v3.*w1.^2+ ...
        (-36).*a21.^2.*a22.*v3.*w1.^2+48.*a22.*d31.*v3.*w1.^2+72.*a21.* ...
        a22.*d31.*v3.*w1.^2+(-36).*a22.*d31.^2.*v3.*w1.^2+(-12).*vv3.* ...
        w1.^2+(-48).*a21.*vv3.*w1.^2+(-36).*a21.^2.*vv3.*w1.^2+48.*d31.* ...
        vv3.*w1.^2+72.*a21.*d31.*vv3.*w1.^2+(-36).*d31.^2.*vv3.*w1.^2+( ...
        -12).*a21.*v3.*w1.^3+(-24).*a21.^2.*v3.*w1.^3+(-12).*a21.^3.*v3.* ...
        w1.^3+12.*d31.*v3.*w1.^3+48.*a21.*d31.*v3.*w1.^3+36.*a21.^2.*d31.* ...
        v3.*w1.^3+(-24).*d31.^2.*v3.*w1.^3+(-36).*a21.*d31.^2.*v3.*w1.^3+ ...
        12.*d31.^3.*v3.*w1.^3+(-1).*b1.*w1.^4);
    aa21=0.0;
    aa22=(1/6).*w1.^(-1).*(a22.^3+3.*a22.^2.*w1+3.*a21.*a22.^2.*w1+(-3).* ...
        a22.^2.*d31.*w1+6.*a21.*a22.*w1.^2+3.*a21.^2.*a22.*w1.^2+(-6).* ...
        a22.*d31.*w1.^2+(-6).*a21.*a22.*d31.*w1.^2+3.*a22.*d31.^2.*w1.^2+ ...
        3.*a21.*w1.^3+3.*a21.^2.*w1.^3+a21.^3.*w1.^3+(-2).*d31.*w1.^3+(-6) ...
        .*a21.*d31.*w1.^3+(-3).*a21.^2.*d31.*w1.^3+3.*d31.^2.*w1.^3+3.* ...
        a21.*d31.^2.*w1.^3+(-1).*d31.^3.*w1.^3);
    
    d=[1   0;
        0    1 ;
        d31  d32;];
    b=[b1,  b2];
    A=[0  0 0;
        0  0 0;
        w1*a21  a22  0;];
    Ahat =[0  0 0;
        0  0 0;
        w1*w1*aa21  aa22  0;];
    v=[w1*v1,  v2, v3] ;
    vhat=[w1*w1*vv1,  vv2,  vv3] ;
    
end
%% TDMSRK2239
if chars==2239
    xx= [0.0930468386721641,0.0815331756382850,0.520345787428202,0.541558149360853,-1.14121445587939];
    b1 = xx(1);   v1 = xx(2); v3 = xx(3);   a22 = xx(4);  w1=w(1);
    a22=0.5416*w1^(-0.04681) ;
    b1=  0.1397*exp(0.8707*w1) -0.2405*w1^1.319;
    v1=  0.6519  *exp(-2.359 *w1) +0.01983*w1^1.807 ;
    v3=   0.3203  *exp(  0.4813 *w1)  *w1^( -0.2711) ;
    %  [0.0930468386721641,0.0815331756382850,0.520345787428202,0.541558149360853,-1.14121445587939]
    b2=1 - b1 ;
    d31 = 0; vv1 = 0; aa21 =  0; a21 = 0; d32=1-d31;
    v2=1+(-1).*v3+b1.*w1+(-1).*v1.*w1;
    
    vv2=(-1/6).*a22.^(-1).*(1+3.*a22.*((-1)+a22.*v3)+3.*a22.*(b1+(-2).*v1) ...
        .*w1.^2+(b1+(-3).*v1).*w1.^3);
    vv3=(1/6).*a22.^(-1).*(1+(-3).*a22.^2.*v3+(b1+(-3).*v1).*w1.^3);
    aa22=(1/2).*a22.^2;
    
    d=[1   0;
        0    1 ;
        d31  d32;];
    b=[b1,  b2];
    A=[0  0 0;
        0  0 0;
        w1*a21  a22  0;];
    Ahat =[0  0 0;
        0  0 0;
        w1*w1*aa21  aa22  0;];
    v=[w1*v1,  v2, v3] ;
    vhat=[w1*w1*vv1,  vv2,  vv3] ;
end
%% TDMSRK224
L=1;
if chars==224

    b1 = x(1);   v3 = x(2); vv3 = x(3); a21 = x(4);  d31 = x(5);  w1=w(1);
    %     b1 = x(1);   v3 = x(2); vv3 = x(3); a21 = x(4);  a22 = x(5); d31 = x(6);  w1=w(1);
    aa21=0;
    a22= (-a21)*w1 + d31*w1 + (3*a21*w1^3 - d31*w1^3)^(1/3);
    
    b2=1 - b1 ;      d32=1 - d31;
    %     [1,0.0776580383519972,0.301371373051056,0.199267993301564,
    %      0.0907759259244949,0.566387706162208,0.0859182382169720,-0.946487048652487]
    %     [0.0776580383519972,0.301371373051056,0.199267993301564,
    %         0.0907759259244949,0.0859182382169720,-0.946487048652487]

    v1=(1/2).*w1.^(-4).*(1+(-4).*a22.^3.*v3+(-12).*a22.^2.*vv3+2.*w1+(-6) ...
        .*a22.^2.*v3.*w1+(-12).*a21.*a22.^2.*v3.*w1+12.*a22.^2.*d31.*v3.* ...
        w1+(-12).*a22.*vv3.*w1+(-24).*a21.*a22.*vv3.*w1+24.*a22.*d31.* ...
        vv3.*w1+(-12).*a21.*a22.*v3.*w1.^2+(-12).*a21.^2.*a22.*v3.*w1.^2+ ...
        12.*a22.*d31.*v3.*w1.^2+24.*a21.*a22.*d31.*v3.*w1.^2+(-12).*a22.* ...
        d31.^2.*v3.*w1.^2+(-12).*a21.*vv3.*w1.^2+(-12).*a21.^2.*vv3.* ...
        w1.^2+12.*d31.*vv3.*w1.^2+24.*a21.*d31.*vv3.*w1.^2+(-12).*d31.^2.* ...
        vv3.*w1.^2+(-6).*a21.^2.*v3.*w1.^3+(-4).*a21.^3.*v3.*w1.^3+12.* ...
        a21.*d31.*v3.*w1.^3+12.*a21.^2.*d31.*v3.*w1.^3+(-6).*d31.^2.*v3.* ...
        w1.^3+(-12).*a21.*d31.^2.*v3.*w1.^3+4.*d31.^3.*v3.*w1.^3+b1.*  w1.^4);
    v2=(1/2).*w1.^(-3).*((-1)+4.*a22.^3.*v3+12.*a22.^2.*vv3+(-2).*w1+6.* ...
        a22.^2.*v3.*w1+12.*a21.*a22.^2.*v3.*w1+(-12).*a22.^2.*d31.*v3.*w1+ ...
        12.*a22.*vv3.*w1+24.*a21.*a22.*vv3.*w1+(-24).*a22.*d31.*vv3.*w1+ ...
        12.*a21.*a22.*v3.*w1.^2+12.*a21.^2.*a22.*v3.*w1.^2+(-12).*a22.* ...
        d31.*v3.*w1.^2+(-24).*a21.*a22.*d31.*v3.*w1.^2+12.*a22.*d31.^2.* ...
        v3.*w1.^2+12.*a21.*vv3.*w1.^2+12.*a21.^2.*vv3.*w1.^2+(-12).*d31.* ...
        vv3.*w1.^2+(-24).*a21.*d31.*vv3.*w1.^2+12.*d31.^2.*vv3.*w1.^2+2.* ...
        w1.^3+(-2).*v3.*w1.^3+6.*a21.^2.*v3.*w1.^3+4.*a21.^3.*v3.*w1.^3+( ...
        -12).*a21.*d31.*v3.*w1.^3+(-12).*a21.^2.*d31.*v3.*w1.^3+6.* ...
        d31.^2.*v3.*w1.^3+12.*a21.*d31.^2.*v3.*w1.^3+(-4).*d31.^3.*v3.* w1.^3+b1.*w1.^4);
    vv1=(1/12).*w1.^(-4).*(3+(-12).*a22.^3.*v3+(-36).*a22.^2.*vv3+4.*w1+( ...
        -12).*a22.^2.*v3.*w1+(-36).*a21.*a22.^2.*v3.*w1+36.*a22.^2.*d31.* ...
        v3.*w1+(-24).*a22.*vv3.*w1+(-72).*a21.*a22.*vv3.*w1+72.*a22.*d31.* ...
        vv3.*w1+(-24).*a21.*a22.*v3.*w1.^2+(-36).*a21.^2.*a22.*v3.*w1.^2+ ...
        24.*a22.*d31.*v3.*w1.^2+72.*a21.*a22.*d31.*v3.*w1.^2+(-36).*a22.* ...
        d31.^2.*v3.*w1.^2+(-24).*a21.*vv3.*w1.^2+(-36).*a21.^2.*vv3.* ...
        w1.^2+24.*d31.*vv3.*w1.^2+72.*a21.*d31.*vv3.*w1.^2+(-36).*d31.^2.* ...
        vv3.*w1.^2+(-12).*a21.^2.*v3.*w1.^3+(-12).*a21.^3.*v3.*w1.^3+24.* ...
        a21.*d31.*v3.*w1.^3+36.*a21.^2.*d31.*v3.*w1.^3+(-12).*d31.^2.*v3.* ...
        w1.^3+(-36).*a21.*d31.^2.*v3.*w1.^3+12.*d31.^3.*v3.*w1.^3+b1.* w1.^4);
    vv2=(1/12).*w1.^(-2).*(3+(-12).*a22.^3.*v3+(-36).*a22.^2.*vv3+8.*w1+( ...
        -24).*a22.^2.*v3.*w1+(-36).*a21.*a22.^2.*v3.*w1+36.*a22.^2.*d31.* ...
        v3.*w1+(-48).*a22.*vv3.*w1+(-72).*a21.*a22.*vv3.*w1+72.*a22.*d31.* ...
        vv3.*w1+6.*w1.^2+(-12).*a22.*v3.*w1.^2+(-48).*a21.*a22.*v3.*w1.^2+ ...
        (-36).*a21.^2.*a22.*v3.*w1.^2+48.*a22.*d31.*v3.*w1.^2+72.*a21.* ...
        a22.*d31.*v3.*w1.^2+(-36).*a22.*d31.^2.*v3.*w1.^2+(-12).*vv3.* ...
        w1.^2+(-48).*a21.*vv3.*w1.^2+(-36).*a21.^2.*vv3.*w1.^2+48.*d31.* ...
        vv3.*w1.^2+72.*a21.*d31.*vv3.*w1.^2+(-36).*d31.^2.*vv3.*w1.^2+( ...
        -12).*a21.*v3.*w1.^3+(-24).*a21.^2.*v3.*w1.^3+(-12).*a21.^3.*v3.* ...
        w1.^3+12.*d31.*v3.*w1.^3+48.*a21.*d31.*v3.*w1.^3+36.*a21.^2.*d31.* ...
        v3.*w1.^3+(-24).*d31.^2.*v3.*w1.^3+(-36).*a21.*d31.^2.*v3.*w1.^3+ ...
        12.*d31.^3.*v3.*w1.^3+(-1).*b1.*w1.^4);
%     aa21=(1/6).*w1.^(-3).*((-1).*a22.^3+(-3).*a21.*a22.^2.*w1+3.*a22.^2.* ...
%         d31.*w1+(-3).*a21.^2.*a22.*w1.^2+6.*a21.*a22.*d31.*w1.^2+(-3).* ...
%         a22.*d31.^2.*w1.^2+3.*a21.*w1.^3+(-1).*a21.^3.*w1.^3+(-1).*d31.* ...
%         w1.^3+3.*a21.^2.*d31.*w1.^3+(-3).*a21.*d31.^2.*w1.^3+d31.^3.*   w1.^3);
    aa22=(1/6).*w1.^(-1).*(a22.^3+3.*a22.^2.*w1+3.*a21.*a22.^2.*w1+(-3).* ...
        a22.^2.*d31.*w1+6.*a21.*a22.*w1.^2+3.*a21.^2.*a22.*w1.^2+(-6).* ...
        a22.*d31.*w1.^2+(-6).*a21.*a22.*d31.*w1.^2+3.*a22.*d31.^2.*w1.^2+ ...
        3.*a21.*w1.^3+3.*a21.^2.*w1.^3+a21.^3.*w1.^3+(-2).*d31.*w1.^3+(-6) ...
        .*a21.*d31.*w1.^3+(-3).*a21.^2.*d31.*w1.^3+3.*d31.^2.*w1.^3+3.* ...
        a21.*d31.^2.*w1.^3+(-1).*d31.^3.*w1.^3);
    
    d=[1   0;
        0    1 ;
        d31  d32;];
    b=[b1,  b2];
    A=[0  0 0;
        0  0 0;
        w1*a21  a22  0;];
    Ahat =[0  0 0;
        0  0 0;
        w1*w1*aa21  aa22  0;];
    v=[w1*v1,  v2, v3] ;
    vhat=[w1*w1*vv1,  vv2,  vv3] ;
    
end
%% TDMSRK223
if chars==223
    b1 = x(1);   v1 = x(2); v3 = x(3);   a22 = x(4);  w1=w(1);
%  [0.0930468386721641,0.0815331756382850,0.520345787428202,0.541558149360853,-1.14121445587939]
    b2=1 - b1 ;
    d31 = 0; vv1 = 0; aa21 =  0; a21 = 0; d32=1-d31;
    v2=1+(-1).*v3+b1.*w1+(-1).*v1.*w1;
    
    vv2=(-1/6).*a22.^(-1).*(1+3.*a22.*((-1)+a22.*v3)+3.*a22.*(b1+(-2).*v1) ...
        .*w1.^2+(b1+(-3).*v1).*w1.^3);
    vv3=(1/6).*a22.^(-1).*(1+(-3).*a22.^2.*v3+(b1+(-3).*v1).*w1.^3);
    aa22=(1/2).*a22.^2;
    
    d=[1   0;
        0    1 ;
        d31  d32;];
    b=[b1,  b2];
    A=[0  0 0;
        0  0 0;
        w1*a21  a22  0;];
    Ahat =[0  0 0;
        0  0 0;
        w1*w1*aa21  aa22  0;];
    v=[w1*v1,  v2, v3] ;
    vhat=[w1*w1*vv1,  vv2,  vv3] ;
end

%% TDMS34
if chars==34
    b1 = x(1);   b2 = x(2);   v2 = x(3);   v3 = x(4);     w1=w(1); w2=w(2);
    b3=1 - b1 - b2;
    
    v1=w1.^(-1).*(1+(-1).*v3+b1.*w1+(b1+b2+(-1).*v2).*w2);
    vv1=(1/12).*w1.^(-3).*(w1+w2).^(-1).*(1+3.*b1.*w1.^4+2.*w2+(-1).* ...
        w2.^3.*(2+(-2).*v3+(b1+b2).*w2)+6.*w1.^2.*w2.*(1+(-1).*v3+(b1+b2+( ...
        -1).*v2).*w2)+(-4).*w1.^3.*((-1)+v3+((-2).*b1+(-1).*b2+v2).*w2));
    vv2=1/12 *w1.^(-1).*w2.^(-3).*((-1)+b1.*w1.^4+(-2).*w2+w2.^3.*(2+( ...
        -2).*v3+(b1+b2).*w2)+6.*w1.^2.*w2.*(1+(-1).*v3+(b1+b2+(-1).*v2).* ...
        w2)+(-2).*w1.^3.*((-1)+v3+((-2).*b1+(-1).*b2+v2).*w2)+2.*w1.*((-1) ...
        +w2.^2.*(3+(-3).*v3+2.*(b1+b2).*w2)));
    vv3=(1/12).*w2.^(-1).*(w1+w2).^(-1).*(1+(-1).*b1.*w1.^4+2.*w1.^3.*(( ...
        -1)+v3+(-1).*(b1+b2+(-1).*v2).*w2)+w2.*(4+w2.*(6+w2.*(4+(-4).*v3+( ...
        b1+b2).*w2)))+2.*w1.*(1+w2.*(3+w2.*(3+(-3).*v3+(b1+b2).*w2))));
    
    d=[1  0 0;
        0 1 0;
        0 0 1;];
    b=[b1,  b2,  b3];
    A=[0  0 0;
        0  0 0;
        0  0 0;];
    Ahat =[0  0 0;
        0  0 0;
        0  0 0;];
    v=[w1*v1, w2*v2, v3] ;
    vhat=[w1*w1*vv1, w2*w2*vv2, vv3] ;
end
%% TDMS341
if chars==341
   % 0.0265839426174917	0.347059703965716	-0.429361803536009
    b1 = x(1);   b2 = x(2);  w1=w(1); w2=w(2);
    b3=1 - b1 - b2;
    v1 = 0; vv2 = 0;

vv1 = (3 + 4*w2 + b2*w2^4 - b1*(3*w1 - w2)*(w1 + w2)^3)/...
       (12*w1^2*(w1 + w2)*(3*w1 + w2)); 
  v3 = (-2 - 3*w2 + w2^3 + 6*vv1*w1^2*(w1 + w2)*...
            (4*w1 + w2) + b1*w1^2*(w1 + w2)*(2*w1 + 3*w2))/w2^3; 
  v2 = (-1 - w2 + b1*w1*(w1 + w2)^3 + 6*vv1*w1^2*(w1 + w2)*  (2*w1 + w2))/w2^4; 
vv3 = (-vv1)*w1^2 + v2*w2^2 + (1/2)*(1 - b1*(-w1 - w2)^2 -    b2*w2^2); 
    
    d=[1  0 0;
        0 1 0;
        0 0 1;];
    b=[b1,  b2,  b3];
    A=[0  0 0;
        0  0 0;
        0  0 0;];
    Ahat =[0  0 0;
        0  0 0;
        0  0 0;];
    v=[w1*v1, w2*v2, v3] ;
    vhat=[w1*w1*vv1, w2*w2*vv2, vv3] ;
end

%% TDMS349
if chars==349
    b1 = x(1);   b2 = x(2);  w1=w(1); w2=w(2);
% b1=a + b*x^(c)*exp(c1*x) *y^(f)
% b2=a + b*x^(c) *y^(f) 

b1= 22.27   + -22.26*w1^( 0.00209 )*exp(-0.0006788*w1) *w2^( -0.000214 );
b2= 0.06399   + 0.2848*w1^(-0.2672) *w2^(  -1.599 ) ;
    b3=1 - b1 - b2;
    v1 = 0; vv2 = 0;

vv1 = (3 + 4*w2 + b2*w2^4 - b1*(3*w1 - w2)*(w1 + w2)^3)/...
       (12*w1^2*(w1 + w2)*(3*w1 + w2)); 
  v3 = (-2 - 3*w2 + w2^3 + 6*vv1*w1^2*(w1 + w2)*...
            (4*w1 + w2) + b1*w1^2*(w1 + w2)*(2*w1 + 3*w2))/w2^3; 
  v2 = (-1 - w2 + b1*w1*(w1 + w2)^3 + 6*vv1*w1^2*(w1 + w2)*  (2*w1 + w2))/w2^4; 
vv3 = (-vv1)*w1^2 + v2*w2^2 + (1/2)*(1 - b1*(-w1 - w2)^2 -    b2*w2^2); 
    
    d=[1  0 0;
        0 1 0;
        0 0 1;];
    b=[b1,  b2,  b3];
    A=[0  0 0;
        0  0 0;
        0  0 0;];
    Ahat =[0  0 0;
        0  0 0;
        0  0 0;];
    v=[w1*v1, w2*v2, v3] ;
    vhat=[w1*w1*vv1, w2*w2*vv2, vv3] ;
end

%% TDMS340
if chars==340
    b1 = x(1);   b2 = x(2);   v1= x(3);  v2 = x(4);   v3 = x(5);
    vv1= x(6);  vv2 = x(7);   vv3 = x(8);  %v1=0;  vv2=0;
    w1=w(1); w2=w(2);
    L=[-w1-w2;
        -w2;
        0 ];
    b3=1 - b1 - b2;
    d=[1  0 0;
        0 1 0;
        0 0 1;];
    b=[b1,  b2,  b3];
    A=[0  0 0;
        0  0 0;
        0  0 0;];
    Ahat =[0  0 0;
        0  0 0;
        0  0 0;];
    v=[w1*v1, w2*v2, v3] ;
    vhat=[w1*w1*vv1, w2*w2*vv2, vv3] ;
end
%% TDMS230
if chars==230
    
    b1 = x(1);  v1 = x(2);  v2 = x(3);   vv1 = x(4);  vv2 = x(5);    w1=w(1);
    b2=1-b1;
    L=[-w1;
        0 ];
    d=[1, 0
        0, 1];
    b=[b1,  b2];
    A=[0  0;
        0  0;];
    Ahat =[0  0;
        0  0;];
    v=[w1*v1, v2] ;
    vhat=[w1*w1*vv1, vv2] ;
end
%% TDMS240
if chars==240
    b1 = x(1);  v1 = x(2);  v2 = x(3);   vv1 = x(4);  vv2 = x(5);    w1=w(1);
    b2=1-b1;
    L=[-w1;
        0 ];
    d=[1, 0
        0, 1];
    b=[b1,  b2];
    A=[0  0;
        0  0;];
    Ahat =[0  0;
        0  0;];
    v=[w1*v1, v2] ;
    vhat=[w1*w1*vv1, vv2] ;
end
%% TDMS23
if chars==23
    b1 = x(1);      w1=w(1);
    b2=1-b1;
     vv1 = 0;
    v1 = 1/3* (b1 + 1/w1^3);
    v2 = 1 - 1/(3 *w1^2) + (2 *b1 *w1)/3; 
    vv2 = ( 2 + 3* w1 - b1* w1^3)/(6* w1);
    d=[1, 0
        0, 1];
    b=[b1,  b2];
    A=[0  0;
        0  0;];
    Ahat =[0  0;
        0  0;];
    v=[w1*v1, v2] ;
    vhat=[w1*w1*vv1, vv2] ;
end
%% TDMS239
if chars==239
    b1 = x(1);      w1=w(1);

    b1= 0.2013 *w1^-1.592 ;
    b2=1-b1;
     vv1 = 0;
    v1 = 1/3* (b1 + 1/w1^3);
    v2 = 1 - 1/(3 *w1^2) + (2 *b1 *w1)/3; 
    vv2 = ( 2 + 3* w1 - b1* w1^3)/(6* w1);
    d=[1, 0
        0, 1];
    b=[b1,  b2];
    A=[0  0;
        0  0;];
    Ahat =[0  0;
        0  0;];
    v=[w1*v1, v2] ;
    vhat=[w1*w1*vv1, vv2] ;
end
%% TDMS231
if chars==231
    
    b1 = x(1);   v2 = x(2);      w1=w(1);
    b2=1-b1;
    v1=b1+(1+(-1).*v2).*w1.^(-1);
    vv1=(1/6).*w1.^(-3).*((-1)+w1.^2.*(3+(-3).*v2+2.*b1.*w1));
    vv2=(1/6).*w1.^(-1).*(1+w1.*(3+w1.*(3+(-3).*v2+b1.*w1)));
    d=[1, 0
        0, 1];
    b=[b1,  b2];
    A=[0  0;
        0  0;];
    Ahat =[0  0;
        0  0;];
    v=[w1*v1, v2] ;
    vhat=[w1*w1*vv1, vv2] ;
end
%% TDMSRK23
% if s==23
% a31 = x(1);    a32 = x(2);     aa31 = x(3);   aa32 = x(4);
% v1 = x(5); v2 = x(6);  v3 = x(7);
% vv1 = x(8); vv2 = x(9);  vv3 = x(10);
% d31 = x(11);  d32 = x(12); b1 = x(13);  b2 = x(14);
% %  a21 = 0.5; v2 = 0.38; ww1 = 0.0; ww2 = 0.0; w1 = 0.0; w2 = 0.1;
%
% d=[1, 0
%    0, 1
%    d31, d32];
% b=[b1,  b2];
%
% A=[0 0 0;
%    0 0 0;
%    a31,a32, 0;];
%
% Ahat=[0 0 0;
%       0 0 0;
%       aa31,aa32, 0;];
% v=[v1,v2,v3] ;
% vhat=[vv1,vv2,vv3] ;
% end
%

end
% count1= (stage-1)*step;  %d
% count2= count1 + step;   %b
% count3= count2 + (2*step+stage-2)*(stage-1)/2; %A
% count4= count3 + (step+stage-1);    %v
% count5 = count4 +(2*step+stage-2)*(stage-1)/2;  %Ahat
% count6= count5 + (step+stage-1);     %vhat
% di=eye(step);
% d=[di;reshape(x(1: count1),stage-1,step)];
% b= reshape(x(count1+1:count2),1,step) ;
% A=zeros(step+stage-1, step+stage-1);
% count=count2;
% for i = step+1 : step+stage-1
% A(i,1:i-1)=x(count+1:count+i-1);
% count= count+i-1;
% end
% v=[ reshape(x(count3+1:count4),1,step+stage-1)];
% Ahat=zeros(step+stage-1, step+stage-1);
% count=count4;
% for i = step+1 : step+stage-1
% Ahat(i,1:i-1)=x(count+1:count+i-1);
% count= count+i-1;
% end
% vhat=[ reshape(x(count5+1:count6),1,step+stage-1)];
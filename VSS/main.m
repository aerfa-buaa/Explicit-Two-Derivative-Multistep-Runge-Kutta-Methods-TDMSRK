clear all
% close all
warning off
clc


K=1/sqrt(2);
K=4
chars=223;  %    选择具体的时间方法
%%
 if chars==223
    wq=1.0;
    for ii =1:5
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.00
    end
 end
 if chars==224
    wq=1.0;
    for ii =1:51
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.00
    end
 end
 
if chars==34
    a_jj=1;
    wq1=1.0;
    for ii =1:1
        wq2=1.0;
        for jj =1:a_jj
            ii
            jj
            w=[wq1, wq2];
            all_coff((ii-1)*a_jj+jj,:)=[w,opt_mdrk(chars, w ,K)];
            wq2=wq2+0.1
        end
        wq1=wq1+0.1
    end
end
 if chars==23
    wq=0.6;
    for ii =1:51
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.02
    end
 end
  if chars==230
    wq=0.6;
    for ii =1:51
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.02
    end
  end
    if chars==240
    wq=1;
    for ii =1:11
        ii
        w=[wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.0
    end
end
 if chars==340
    wq=1.0;
    for ii =1:41
        ii
        w=[wq,wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.01 
    end
end

 if chars==341
    wq=1.0;
    for ii =1:41
        ii
        w=[wq,wq];
        all_coff(ii,:)=[w,opt_mdrk(chars, w ,K)];
        wq=wq+0.0  
    end
end


% k1=[0.1, 0.2, 0.3,  0.4,  0.5, 0.6, 1/sqrt(2), 0.8  1,  1.5,  2,  3,  4]; ii=0;
% [hh,nb0]=size(k1);
% nb1=40;
%
% method=25;
% for iqq0=1:nb1
% %           K=k1(iqq0);
%           K=0.1*iqq0
%        opt_r =opt_mdrk(method,K);
%         ii=ii+1;
%         all_coff(ii,1)=K;
%         [hh21,r32]=size(opt_r);
%         for  ci=1:r32
%         all_coff(ii,1+ci)= opt_r(ci) ;
%         end
%     end

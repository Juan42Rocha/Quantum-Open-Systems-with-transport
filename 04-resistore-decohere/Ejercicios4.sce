L=100;
E=0.0;
nu=0.0;
OnSiteE=0.0;

transfile=mopen('Resistivity_nu:'+string(nu)+'_ContactosFicticios.ou','wt');


//For Anderson localization:

  rho=zeros(1,1);
  aux=1;

  for LL=5:5:L
    auxM= zeros(LL,LL);
    TDAlpha=0.0; TBetaS=0.0;
    T=zeros(LL,LL); 
    TT=0.0; TTCont=0.0;
    for j=1:50
      H= diag(ones(LL-1,1),1).*.1 + diag(ones(LL-1,1),-1).*.1 + diag(grand(LL,1,'nor',0,1));

      eta=  1*(E/2.0-%i*sqrt(1-E**2/4.0));
      
      sigCont= eye(LL,LL)*(-%i*nu);
      
      sigS=zeros(LL,LL);
      sigD=sigS;
      sigS(1,1)=eta;
      sigD(LL,LL)=eta;
      EM=eye(LL,LL)*E;

      G=(EM-H-sigCont-sigS-sigD)**(-1);
      
      gammaD= %i*(sigD-sigD');
      gammaS= %i*(sigS-sigS');

      for alpha=1:LL
        auxM2= auxM; auxM2(alpha,alpha)= (-%i*nu);
      
        TDAlpha(alpha)= sum( diag( gammaD* (G*  (%i*(auxM2-auxM2')*G') ) ));
        TBetaS(alpha)= sum( diag(%i*(auxM2-auxM2')* (G* (gammaS*G'))  ));
        for beta=1:LL
          auxM3= auxM; auxM3(beta,beta)= (-%i*nu);
      
          if alpha == beta then
            for k=1:LL
              auxM4= auxM; auxM4(k,k)= (-%i*nu);
              if alpha ~= k then
                T(alpha,beta)= T(alpha,beta) + sum( diag(%i*(auxM2-auxM2')* (G* (%i*(auxM4-auxM4')*G'))  ));
              end
            end
            T(alpha,beta)= T(alpha,beta) + sum( diag(%i*(auxM2-auxM2')* (G* (gammaS*G'))  )) + sum( diag(%i*(auxM2-auxM2')* (G* (gammaD*G'))  ));
          else
            T(alpha,beta)= - sum( diag(%i*(auxM2-auxM2')* (G* (%i*(auxM3-auxM3')*G'))  ));
          end
        end
      end
      
      Tao=T**(-1);
      TTCont= TDAlpha'*(Tao*TBetaS);
      
      TT= TT + 4.0*sum(diag(imag(sigS)*G*imag(sigD)*G')) + TTCont;
    end
    TT= TT/50.0;
    Ttotal(aux)=TT;
    aux=aux+1;
  end
 
 

  for i=1:length(Ttotal)
    mfprintf(transfile,'%g   %g\n',5+5*i,1.0/Ttotal(i));
  end
  mfprintf(transfile,'\n\n\n');

mclose('all');
exit

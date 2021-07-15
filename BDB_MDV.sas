/*MDV: minimal detective value*/

%macro simmdv(n=, eventnumber=, w=, miu_Iprior= ,tau_Iprior=0.173,miu_Vprior=0,tau_Vprior=2.59);

     data postmdv;
     do HR = 1 to 0.6  by -0.001; 
        datamean=log(HR);
        dataSE2= 4/&eventnumber.;

        /*mean and sd of  difference from posterior based on informative prior ; */
        rou_I=1/&tau_Iprior**2+1/dataSE2;
        miu_post_I=(1/rou_I)*(&miu_Iprior/&tau_Iprior**2+datamean/dataSE2);
        sd_post_I=sqrt(1/rou_I);
        /*marginal distribution based on  informative prior*/
        f_I=exp(-0.5*(datamean-&miu_Iprior)**2/(dataSE2+&tau_Iprior**2))/sqrt(dataSE2+&tau_Iprior**2);

        /*mean and sd of  difference from posterior based on vague prior ; */
        rou_V=1/&tau_Vprior**2+1/dataSE2;
        miu_post_V=(1/rou_V)*(&miu_Vprior/&tau_Vprior**2+datamean/dataSE2);
        sd_post_V=sqrt(1/rou_V);
        /*marginal distribution based on vague prior*/
        f_V=exp(-0.5*(datamean-&miu_Vprior)**2/(dataSE2+&tau_Vprior**2))/sqrt(dataSE2+&tau_Vprior**2);

        /*update weight*/
        wei_post=(&w.*f_I)/(&w.*f_I+(1-&w.)*f_V);

        Prob=wei_post*CDF("normal",0 ,miu_post_I,sd_post_I) +(1-wei_post)*CDF("normal",0,miu_post_V,sd_post_V);
        /*The MDV is the max HR with Prob>xx*/
        if Prob>0.95 then success=1;
        else success=0;
        output;
     end;
  run;

  data mdv;
    set postmdv(where=(success=1) keep=hr success);
    if _n_=1;
    call symput ('MDV',hr);
    n=&n.;
    eventnumber=&eventnumber.;
    miu_Iprior=&miu_Iprior.;
    tau_Iprior=&tau_Iprior.;
    miu_Vprior=&miu_Vprior.;
    tau_Vprior=&tau_Vprior.;
  run;
%put %str(MDV=)&mdv. %str(, prior=)&miu_Iprior;
%mend;



%simmdv(n=100, eventrate=0.5, w=0.4, miu_Iprior=log(0.5),tau_Iprior=0.173,miu_Vprior=0,tau_Vprior=2.59);




/*
* Purpose:
* Input:
* Output:
* 
*/
dm log 'clear;output;';
options source notes nocenter nonumber nodate mprint;
%let Ntrial=50000;

%macro mydata(eventnumber, diff_start,diff_end,diff_step);
data mydata;
  do w=0.35 to 0.5 by 0.05;
  do n=80 to 112 by 4;
    do diff=&diff_start to &diff_end by &diff_step;
      do trial=1 to &Ntrial;
          /*sample the mean using calculated vague prior SD*/
          dataSE2= 4/&eventnumber. ;
          datamean=rand("NORMAL", log(diff), sqrt(dataSE2));
          output;
       end;
     end;
    do diff=1;
       do trial=1 to &Ntrial;
          /*sample the mean using calculated vague prior SD*/
          dataSE2= 4/&eventnumber. ;
          datamean=rand("NORMAL", log(diff), sqrt(dataSE2));
          output;
       end;
    end;
   end;
   end;
run;
%mend;

%mydata(0.5, diff_start=0.5,diff_end=0.6,diff_step=0.1);


%macro sim(set, miu_Iprior,tau_Iprior,miu_Vprior,tau_Vprior,n_start,n_end,n_step,wei_start,wei_end,wei_step);
*I: informative prior  V: vague prior;
data posterior;
  set &set;
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
        /*method1 accordting to Schmidli et al literature*/
        wei_post=(w*f_I)/(w*f_I+(1-w)*f_V);
                
        /*for non-inferiority,just change the 0 to -0.1 */
        Prob1=wei_post*CDF("normal",0 ,miu_post_I,sd_post_I) +(1-wei_post)*CDF("normal",0,miu_post_V,sd_post_V);
        if Prob1>0.95 then success=1; else success=0;

  run;
%mend;

*%sim(set=mydata,miu_Iprior=&prior_mean,tau_Iprior=&prior_std,miu_Vprior=0,tau_Vprior=2?);
%sim(set=mydata,miu_Iprior=log(0.6),tau_Iprior=0.173 ,miu_Vprior=0,tau_Vprior=2);
*%sim(set=mydata,miu_Iprior=log(0.66),tau_Iprior=0.173 ,miu_Vprior=0,tau_Vprior=2.59);
*%sim(set=mydata,miu_Iprior=log(0.72),tau_Iprior=0.173 ,miu_Vprior=0,tau_Vprior=2.59);

ods listing close;
proc means data=posterior mean ;
  var success wei_post;
  class diff w n;
  ods output summary=des;
run;

ods excel file="result.xlsx";

proc print data=des;
run;

ods excel close;


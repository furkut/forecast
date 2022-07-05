import java.util.ArrayList;
import java.util.List;

public class falgo {double[]x;
    public double [] exponentialSmoothing(double[]x){

        double[] exSmoothing = new double[x.length];
        exSmoothing[0] =x[0];
        double a =0.2;
        int i ;
        for (i =1; i < exSmoothing.length; i++) {
            exSmoothing[i] = a*x[i-1]+(1-a)*exSmoothing[i-1];
        }

        return exSmoothing;
    }
    public double[] doubleExponentialSmoothing(double[]x){
        double[] doubleExSmoothing = new double[x.length];
        double[] st = new double[x.length];
        double[] gt = new double[x.length];
        double a = 0.2;
        double b = 0.2;
        st[0] = 200;
        gt[0]= 50;
        int i;
        for (i =1; i < st.length; i++) {
                st[i] = a * x[i-1] + (1 - a) * (st[i - 1] + gt[i - 1]);
                gt[i] = b * (st[i] - st[i - 1]) + (1 - b) * gt[i - 1];
        }
        for (i =0; i < x.length; i++) {
            doubleExSmoothing[i] =st[i]+(i+x.length+1)*gt[i];
        }

        return doubleExSmoothing;
    }
    public double[] regressionAnalysis(double[]x){
        double[] regressionAn = new double[x.length];
        double sumMonth_Sales=0;
        double sumSales=0;
        int i;
        double monthsSum =0;
        int monthSumSqrt=0;
        for (i =0; i < x.length; i++) {
            sumSales +=x[i];
            sumMonth_Sales+=(i+1)*x[i];
            monthsSum +=i+1;
            monthSumSqrt+=(i+1)*(i+1);
        }
        double b = ((x.length)*sumMonth_Sales-monthsSum*sumSales)/((x.length)*monthSumSqrt-monthsSum*monthsSum);
        double a =sumSales/(x.length)-b*(monthsSum/(x.length));
        for (i =0; i < x.length; i++) {
            regressionAn[i]=a+b*(i+x.length+1);
        }

        return regressionAn;
    }
    public double[] deseasonalizedRegressionAnalysis (double[]x){
        int i;
        double overallavg;
        double sum=0;
        double monthsSum =0;
        double monthSumSqrt=0;
        double deDmdSum=0;
        double deDmd_Period=0;
        double[] deseasonalizedRegAn = new double[x.length];
        double[] pfac = new double[x.length];
        double[] deseasondmd = new double[x.length];

        for (i =0; i < x.length; i++) {
            sum += x[i];
        }

        overallavg=sum/x.length;

        for (i=0;i<x.length/2;i++){
            pfac[i]=((x[i]+x[i+x.length/2])/2)/overallavg;
            deseasondmd[i]=x[i]/pfac[i];
            deseasondmd[i+x.length/2]=x[i+x.length/2]/pfac[i];
        }

        for (i =0; i < x.length; i++) {
            deDmdSum += deseasondmd[i];
            deDmd_Period +=(i+1)*deseasondmd[i];
            monthsSum +=i+1;
            monthSumSqrt+=(i+1)*(i+1);
        }

        double q = monthsSum/x.length;
        double b = ((x.length)*deDmd_Period-monthsSum*deDmdSum)/((x.length)*monthSumSqrt-monthsSum*monthsSum);
        double a=overallavg-b*q;
        for (i =0; i < x.length; i++) {
            deseasonalizedRegAn[i]=a+b*(i+x.length+1);
        }

        return deseasonalizedRegAn;

    }
    public double MSE (double[] exSmoothing,double[] doubleExSmoothing,double[] regressionAn,double[] deseasonalizedRegAn){
        double sumMseEx = 0;
        double sumMseDEx = 0;
        double sumMseReg = 0;
        double sumMseDeReg = 0;
        int z=0;
        int i ;

        for (i =0; i < exSmoothing.length; i++) {
            sumMseEx += x[i]*x[i]-exSmoothing[i]*exSmoothing[i];
            sumMseDEx += x[i]*x[i]-doubleExSmoothing[i]*doubleExSmoothing[i];
            sumMseReg += x[i]*x[i]-regressionAn[i]*regressionAn[i];
            sumMseDeReg += x[i]*x[i]-deseasonalizedRegAn[i]*deseasonalizedRegAn[i];
        }

        double mseEx = sumMseEx/x.length;
        double mseDEx = sumMseDEx/x.length;
        double mseReg = sumMseReg/x.length;
        double mseDeReg = sumMseDeReg/x.length;
        double[] MSE= {mseEx,mseDEx,mseReg,mseDeReg};
        double min = MSE[0] ;

        for (i =0; i < MSE.length; i++) {
            if (min>MSE[i]){
                min = MSE[i];
            }
        }
        for (i =0; i < MSE.length; i++) {
            if (min==MSE[i]){
                z = i;
            }
        }

        return z;
    }

    public double[] min_forecasts(double[] exSmoothing,double[] doubleExSmoothing,double[] regressionAn,double[] deseasonalizedRegAn){
        int i;
        double min=0;
        double[] minNumbers = new double[4];

        for (i =0; i < x.length; i++) {
            if (min>exSmoothing[i]){
                minNumbers[0] = exSmoothing[i];
            }
            if (min>doubleExSmoothing[i]){
                minNumbers[1] = doubleExSmoothing[i];
            }
            if (min>exSmoothing[i]){
                minNumbers[2] = regressionAn[i];
            }
            if (min>exSmoothing[i]){
                minNumbers[3] = deseasonalizedRegAn[i];
            }
        }
        return minNumbers;
    }
    public double[] max_forecasts(double[] exSmoothing,double[] doubleExSmoothing,double[] regressionAn,double[] deseasonalizedRegAn){
        int i;
        double max=0;
        double[] maxNumbers = new double[4];

        for (i =0; i < x.length; i++) {
            if (max<exSmoothing[i]){
                maxNumbers[0] = exSmoothing[i];
            }
            if (max<exSmoothing[i]){
                maxNumbers[1] = doubleExSmoothing[i];
            }
            if (max<exSmoothing[i]){
                maxNumbers[2] = regressionAn[i];
            }
            if (max<exSmoothing[i]){
                maxNumbers[3] = deseasonalizedRegAn[i];
            }
        }
        return maxNumbers;
    }

}

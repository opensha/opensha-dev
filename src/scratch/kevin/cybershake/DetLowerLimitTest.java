package scratch.kevin.cybershake;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.gui.plot.GraphWindow;
import org.opensha.sha.calc.mcer.ASCEDetLowerLimitCalc;
import org.opensha.sha.cybershake.calc.mcer.MCERDataProductsCalc;

public class DetLowerLimitTest {

	public static void main(String[] args) {
		double fa = 1.0; // for site class C
		double fv = 1.3; // for site class C
		double tl = 8;
		
		ArbitrarilyDiscretizedFunc xValsFunc = new ArbitrarilyDiscretizedFunc();
		for (double x=0; x<12d; x+=0.05)
			xValsFunc.set(x, 0d);
		
		DiscretizedFunc lowerLimit = ASCEDetLowerLimitCalc.calc(xValsFunc, fv, fa, tl);
		
		GraphWindow gw = new GraphWindow(lowerLimit, "Deterministic Lower Limit");
		gw.setX_AxisLabel("T (sec)");
		gw.setY_AxisLabel("Sa (g)");
		gw.setDefaultCloseOperation(GraphWindow.EXIT_ON_CLOSE);
	}

}

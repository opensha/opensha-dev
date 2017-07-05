package scratch.peter.ucerf3.calc;

import java.io.File;

import org.opensha.nshmp2.util.Period;

/**
 * Add comments here
 *
 * 
 * @author Peter Powers
 * @version $Id:$
 */
public class AggregateResults {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		File dir = new File(args[0]);
//		Period period = Period.valueOf(Period.class, args[1]);
		
		File dir = new File("/Users/pmpowers/projects/OpenSHA/tmp/UC33/maps/src/UC33brAvg5x_fm31_nobg/");
		Period period = Period.GM0P00;
		UC3_CalcMPJ_MapCompound.aggregateResults(dir, period);
	}

}

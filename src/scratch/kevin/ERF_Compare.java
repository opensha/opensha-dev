package scratch.kevin;

import org.opensha.sha.cybershake.db.MeanUCERF2_ToDB;
import org.opensha.sha.earthquake.ERF;
import org.opensha.sha.earthquake.ProbEqkRupture;
import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;

public class ERF_Compare {
	
	private ERF erf1;
	private ERF erf2;
	
	private boolean probsEqual = true;
	private boolean magsEqual = true;
	private boolean namesEqual = true;
	private boolean rupCountEqual = true;
	
	public ERF_Compare(ERF erf1, ERF erf2) {
		this.erf1 = erf1;
		this.erf2 = erf2;
	}
	
	public void compare() {
		for (int i=0; i<erf1.getNumSources(); i++) {
			ProbEqkSource source1 = erf1.getSource(i);
			ProbEqkSource source2 = erf2.getSource(i);
			
			String name1 = source1.getName();
			String name2 = source2.getName();
			
			if (!name1.equals(name2)) {
				System.out.println("source " + i + " names don't match!");
				System.out.println("'" + name1 + "' vs '" + name2 + "'");
				namesEqual = false;
			}
			
			if (source1.getNumRuptures() != source2.getNumRuptures()) {
				System.out.println("source " + i + " rupts don't match!");
				rupCountEqual = false;
			} else {
				for (int j=0; j<source1.getNumRuptures(); j++) {
					ProbEqkRupture rup1 = source1.getRupture(j);
					ProbEqkRupture rup2 = source2.getRupture(j);
					
					if (rup1.getMag() != rup2.getMag()) {
						magsEqual = false;
					}
					if (rup1.getProbability() != rup2.getProbability()) {
						probsEqual = false;
					}
				}
			}
		}
		
		System.out.println("Names eqal? " + namesEqual);
		System.out.println("Rup count eqal? " + rupCountEqual);
		System.out.println("Mags eqal? " + magsEqual);
		System.out.println("Probs eqal? " + probsEqual);
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ERF erf1 = MeanUCERF2_ToDB.createUCERF2ERF();
		erf1.updateForecast();
		ERF erf2 = MeanUCERF2_ToDB.createUCERF2ERF();
		erf2.getAdjustableParameterList().getParameter(UCERF2.PROB_MODEL_PARAM_NAME)
					.setValue(MeanUCERF2.PROB_MODEL_WGCEP_PREF_BLEND);
		erf2.getTimeSpan().setStartTime(2010);
		erf2.updateForecast();
		
		ERF_Compare comp = new ERF_Compare(erf1, erf2);
		comp.compare();
	}

}

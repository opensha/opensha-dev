package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.opensha.commons.data.Named;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

/**
 * Interface for providing RotD50 results from simulation for analysis and hazard calculation. The rupture type is generic
 * to account for multiple sources (could be a BBP source, RSQSim rupture, CyberShake rupture, etc). Each rupture can have
 * multiple simulations (e.g. a CyberShake rupture variation), which default to equal weight (unless you override 
 * getIndividualSimulationRate(E, double, int, int)).
 * @author kevin
 *
 * @param <E>
 */
public interface SimulationRotDProvider<E> extends Named {
	
	/**
	 * @param site
	 * @param rupture
	 * @param index
	 * @return RotD50 spectrum for the index-th simulation of rupture E at the given site
	 * @throws IOException
	 */
	public DiscretizedFunc getRotD50(Site site, E rupture, int index) throws IOException;
	
	/**
	 * @param site
	 * @param rupture
	 * @return RotD50 spectra for each simulation of rupture E at the given site
	 * @throws IOException
	 */
	public default List<DiscretizedFunc> getRotD50s(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotD50(site, rupture, i));
		return list;
	}
	
	/**
	 * @param site
	 * @param rupture
	 * @param index
	 * @return RotD100 spectrum for the index-th simulation of rupture E at the given site
	 * @throws IOException
	 */
	public DiscretizedFunc getRotD100(Site site, E rupture, int index) throws IOException;
	
	/**
	 * @param site
	 * @param rupture
	 * @return RotD100 spectra for each simulation of rupture E at the given site
	 * @throws IOException
	 */
	public default List<DiscretizedFunc> getRotD100s(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotD100(site, rupture, i));
		return list;
	}
	
	/**
	 * @param site
	 * @param rupture
	 * @return array of [RotD50, RotD100] for the index-th simulation of rupture E at the given site
	 */
	public DiscretizedFunc[] getRotD(Site site, E rupture, int index) throws IOException;

	/**
	 * @param site
	 * @param rupture
	 * @return array of [RotD50, RotD100] for each simulation of rupture E at the given site
	 * @throws IOException
	 */
	public default List<DiscretizedFunc[]> getRotDs(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc[]> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotD(site, rupture, i));
		return list;
	}
	
	/**
	 * @param site
	 * @param rupture
	 * @return ratio of RotD100 to RotD50 for the index-th simulation of rupture E at the given site
	 */
	public DiscretizedFunc getRotDRatio(Site site, E rupture, int index) throws IOException;
	
	/**
	 * @param site
	 * @param rupture
	 * @return ratio of RotD100 to RotD50 for each simulation of rupture E at the given site
	 */
	public default List<DiscretizedFunc> getRotDRatios(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotDRatio(site, rupture, i));
		return list;
	}
	
	/**
	 * @param site
	 * @param rupture
	 * @return the number of simulations available of rupture E at the given site
	 */
	public int getNumSimulations(Site site, E rupture);
	
	/**
	 * @param rupture
	 * @param index
	 * @return the hypocenter of the index-th simulation of rupture E
	 */
	public Location getHypocenter(E rupture, int index);
	
	/**
	 * @param site
	 * @return all ruptures which were simulated at the given site
	 */
	public Collection<E> getRupturesForSite(Site site);
	
	/**
	 * @param site
	 * @return total number of simulations for the given site
	 */
	public default int getNumSimlationsForSite(Site site) {
		int count = 0;
		for (E rupture : getRupturesForSite(site))
			count += getNumSimulations(site, rupture);
		return count;
	}
	
	/**
	 * @return true if RotD50 is available from this source
	 */
	public boolean hasRotD50();
	
	/**
	 * @return true if RotD100 is available from this source
	 */
	public boolean hasRotD100();
	
	/**
	 * Utility method to compute RotD100/RotD50 ratio for an array of [RotD50, RotD100] 
	 * @param spectra
	 * @return
	 */
	public static DiscretizedFunc calcRotDRatio(DiscretizedFunc[] spectra) {
		if (spectra.length == 1)
			// already reducded to a ratio
			return spectra[0];
		Preconditions.checkState(spectra.length == 2);
		ArbitrarilyDiscretizedFunc ret = new ArbitrarilyDiscretizedFunc();
		
		Preconditions.checkState(spectra[0].size() == spectra[1].size());
		
		for (int i=0; i<spectra[0].size(); i++) {
			double period = spectra[0].getX(i);
			Preconditions.checkState((float)period == (float)spectra[1].getX(i));
			double rd50 = spectra[0].getY(i);
			double rd100 = spectra[1].getY(i);
			double ratio = rd100/rd50;
//			System.out.println(period+"s: "+rd100+" / "+rd50+" = "+ratio);
			Preconditions.checkState(Doubles.isFinite(ratio), "Bad ratio! p=%s, rd50=%s, rd100=%s, ratio=%s", period, rd50, rd100, ratio);
			Preconditions.checkState(ratio >= 1d);
			ret.set(period, ratio);
		}
		
		return ret;
	}
	
	/**
	 * @param rupture
	 * @return total annualized rate of the given rupture
	 */
	public double getAnnualRate(E rupture);
	
	/**
	 * Default implementation weights all simulations evenly and returns rupAnnualRate/numSimulations. Can be ovirredden for unequal
	 * simulation weighting.
	 * @param rupture
	 * @param rupAnnualRate
	 * @param simulationIndex
	 * @param numSimulations
	 * @return individual simulation rate of the simulationIndex-th simulation of rupture E
	 */
	public default double getIndividualSimulationRate(E rupture, double rupAnnualRate, int simulationIndex, int numSimulations) {
		return rupAnnualRate/(double)numSimulations;
	}
	
	/**
	 * @param site
	 * @return minimum rate below which curves should be truncated
	 */
	public double getMinimumCurvePlotRate(Site site);
	
	/**
	 * @param rupture
	 * @return magnitude of rupture E
	 */
	public double getMagnitude(E rupture);

}

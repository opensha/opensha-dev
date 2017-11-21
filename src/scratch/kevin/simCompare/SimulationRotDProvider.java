package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.opensha.commons.data.Named;
import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.data.function.DiscretizedFunc;

import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;

public interface SimulationRotDProvider<E> extends Named {
	
	public DiscretizedFunc getRotD50(Site site, E rupture, int index) throws IOException;
	
	public default List<DiscretizedFunc> getRotD50s(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotD50(site, rupture, i));
		return list;
	}
	
	public DiscretizedFunc getRotD100(Site site, E rupture, int index) throws IOException;
	
	public default List<DiscretizedFunc> getRotD100s(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotD100(site, rupture, i));
		return list;
	}
	
	/**
	 * @param site
	 * @param rupture
	 * @return array of [RotD50, RotD100]
	 */
	public DiscretizedFunc[] getRotD(Site site, E rupture, int index) throws IOException;
	
	public default List<DiscretizedFunc[]> getRotDs(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc[]> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotD(site, rupture, i));
		return list;
	}
	
	public DiscretizedFunc getRotDRatio(Site site, E rupture, int index) throws IOException;
	
	public default List<DiscretizedFunc> getRotDRatios(Site site, E rupture) throws IOException {
		ArrayList<DiscretizedFunc> list = new ArrayList<>();
		for (int i=0; i<getNumSimulations(site, rupture); i++)
			list.add(getRotDRatio(site, rupture, i));
		return list;
	}
	
	public int getNumSimulations(Site site, E rupture);
	
	public Collection<E> getRupturesForSite(Site site);
	
	public default int getNumSimlationsForSite(Site site) {
		int count = 0;
		for (E rupture : getRupturesForSite(site))
			count += getNumSimulations(site, rupture);
		return count;
	}
	
	public boolean hasRotD50();
	
	public boolean hasRotD100();
	
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
	
	public double getAnnualRate(E rupture);
	
	public double getMagnitude(E rupture);

}

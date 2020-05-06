package scratch.kevin.simCompare;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.opensha.commons.data.Site;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.imr.param.IntensityMeasureParams.DurationTimeInterval;

import com.google.common.base.Preconditions;

public class SubSetSimulationRotDProvider<E> implements SimulationRotDProvider<E> {
	
	private SimulationRotDProvider<E> simProv;
	private Set<E> ruptures;
	private double rateScalar;
	
	public SubSetSimulationRotDProvider(SimulationRotDProvider<E> simProv, Set<E> ruptures, double rateScalar) {
		this.simProv = simProv;
		this.ruptures = ruptures;
		this.rateScalar = rateScalar;
	}

	@Override
	public String getName() {
		return simProv.getName();
	}

	@Override
	public DiscretizedFunc getRotD50(Site site, E rupture, int index) throws IOException {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getRotD50(site, rupture, index);
	}

	@Override
	public DiscretizedFunc getRotD100(Site site, E rupture, int index) throws IOException {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getRotD100(site, rupture, index);
	}

	@Override
	public DiscretizedFunc[] getRotD(Site site, E rupture, int index) throws IOException {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getRotD(site, rupture, index);
	}

	@Override
	public DiscretizedFunc getRotDRatio(Site site, E rupture, int index) throws IOException {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getRotDRatio(site, rupture, index);
	}

	@Override
	public double getPGV(Site site, E rupture, int index) throws IOException {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getPGV(site, rupture, index);
	}

	@Override
	public double getDuration(Site site, E rupture, DurationTimeInterval interval, int index) throws IOException {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getDuration(site, rupture, interval, index);
	}

	@Override
	public int getNumSimulations(Site site, E rupture) {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getNumSimulations(site, rupture);
	}

	@Override
	public Collection<E> getRupturesForSite(Site site) {
		List<E> siteRups = new ArrayList<>(simProv.getRupturesForSite(site));
		siteRups.retainAll(ruptures);
		return siteRups;
	}

	@Override
	public boolean hasRotD50() {
		return simProv.hasRotD50();
	}

	@Override
	public boolean hasRotD100() {
		return simProv.hasRotD100();
	}

	@Override
	public boolean hasPGV() {
		return simProv.hasPGV();
	}

	@Override
	public boolean hasDurations() {
		return simProv.hasDurations();
	}

	@Override
	public double getAnnualRate(E rupture) {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getAnnualRate(rupture)*rateScalar;
	}

	@Override
	public double getMinimumCurvePlotRate(Site site) {
		return simProv.getMinimumCurvePlotRate(site)*rateScalar;
	}

	@Override
	public double getMagnitude(E rupture) {
		Preconditions.checkState(ruptures.contains(rupture));
		return simProv.getMagnitude(rupture);
	}

	@Override
	public Location getHypocenter(E rupture, int index) {
		return simProv.getHypocenter(rupture, index);
	}

}

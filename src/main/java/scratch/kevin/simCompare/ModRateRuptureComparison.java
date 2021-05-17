package scratch.kevin.simCompare;

import java.util.Collection;

import org.opensha.commons.data.Site;
import org.opensha.sha.earthquake.EqkRupture;

public class ModRateRuptureComparison<E> extends RuptureComparison<E> {

	private RuptureComparison<E> comp;
	private double rate;

	public ModRateRuptureComparison(RuptureComparison<E> comp, double rate) {
		super(comp.getRupture());
		this.comp = comp;
		this.rate = rate;
	}

	@Override
	public Collection<Site> getApplicableSites() {
		return comp.getApplicableSites();
	}

	@Override
	public double getLogMean(Site site, IMT imt) {
		return comp.getLogMean(site, imt);
	}

	@Override
	public double getStdDev(Site site, IMT imt) {
		return comp.getStdDev(site, imt);
	}

	@Override
	public double getDistanceRup(Site site) {
		return comp.getDistanceRup(site);
	}

	@Override
	public double getDistanceJB(Site site) {
		return comp.getDistanceJB(site);
	}

	@Override
	public boolean hasSite(Site site) {
		return comp.hasSite(site);
	}

	@Override
	public boolean isComputed(Site site, IMT imt) {
		return comp.isComputed(site, imt);
	}

	@Override
	public EqkRupture getGMPERupture() {
		return comp.getGMPERupture();
	}

	@Override
	public double getMagnitude() {
		return comp.getMagnitude();
	}

	@Override
	public double getAnnualRate() {
		return rate;
	}

	@Override
	public double getRuptureTimeYears() {
		return comp.getRuptureTimeYears();
	}

}

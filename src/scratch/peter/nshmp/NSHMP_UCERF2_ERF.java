package scratch.peter.nshmp;

import java.util.ArrayList;

import org.opensha.sha.earthquake.ProbEqkSource;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.UCERF2;
import org.opensha.sha.earthquake.rupForecastImpl.WGCEP_UCERF_2_Final.MeanUCERF2.MeanUCERF2;

import com.google.common.collect.Lists;

/**
 * Custom ERF that facilitates comparisons with UCERF3 by not including
 * non-CA b-faults and by clipping grid sources at the border such that ratio
 * maps do not have strong halos around the state border.
 *
 * @author Peter Powers
 */
public class NSHMP_UCERF2_ERF extends MeanUCERF2 {

	private NSHMP08_ClippedGridSourceGen gridSrcGen;
	
	public NSHMP_UCERF2_ERF() {
		gridSrcGen = new NSHMP08_ClippedGridSourceGen();
	}

	@Override
	public ProbEqkSource getSource(int idx) {
		return (idx < allSources.size()) ? 
			allSources.get(idx) : 
			gridSrcGen.getSource(idx - allSources.size());
	}

	@Override
	public int getNumSources() {
		if(backSeisParam.getValue().equals(UCERF2.BACK_SEIS_INCLUDE) ||
				backSeisParam.getValue().equals(UCERF2.BACK_SEIS_ONLY))
			return allSources.size() + gridSrcGen.getNumSources();
		return allSources.size();
	}

	@Override
	protected void mkNonCA_B_FaultSources() {
		// skip non-CA b-faults by returning empty list
		nonCA_bFaultSources = Lists.newArrayList();
	}
	
	@Override
	protected void updateGridSources() {
		// skip adding C-zones to allSurces as they will be handled by
		// gridSrcGen, but they will also be turned off when background
		// siesmicity is turned off
		gridSrcGen.setForecastDuration(timeSpan.getDuration());
	}

	@Override
	public ArrayList<ProbEqkSource> getSourceList(){
		throw new UnsupportedOperationException();
	}

}

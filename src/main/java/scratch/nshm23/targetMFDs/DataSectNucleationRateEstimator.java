package scratch.nshm23.targetMFDs;

import org.opensha.sha.faultSurface.FaultSection;
import org.opensha.sha.magdist.IncrementalMagFreqDist;

public interface DataSectNucleationRateEstimator {
	
	public boolean appliesTo(FaultSection sect);
	
	public double estimateNuclRate(FaultSection sect, IncrementalMagFreqDist curSectSupraSeisMFD);

}

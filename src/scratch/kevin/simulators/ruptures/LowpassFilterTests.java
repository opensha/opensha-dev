package scratch.kevin.simulators.ruptures;

import scratch.kevin.bbp.BBP_Wrapper;
import scratch.kevin.bbp.SeismogramPlotter;
import scratch.kevin.bbp.SpectraPlotter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_Site;

public class LowpassFilterTests {

	public static void main(String[] args) throws IOException {
		File srcFile = new File("/home/kevin/Simulators/catalogs/bruce/rundir2829/event_srfs/event_5304.src");
		File srfFile = new File("/home/kevin/Simulators/catalogs/bruce/rundir2829/event_srfs/event_5304_0.05s_ADJ_VEL.srf");
		
		File outputDir = new File("/tmp/seis_filter");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		List<BBP_Site> sites = new ArrayList<>();
		
		VelocityModel vm = VelocityModel.LA_BASIN_863;
		
		sites.add(new BBP_Site("USC", new Location(34.0192, -118.286), vm.getVs30(),
				RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ));
		
		File sitesFile = new File(outputDir, "sites.stl");
		BBP_Site.writeToFile(sitesFile, sites);
		
		double[] lowpass_freqs = { 00d, 1d, 2d, 3d, 4d, 5d };
		
		for (double lowpassFreq : lowpass_freqs) {
			File freqDir;
			if (lowpassFreq == 0d)
				freqDir = new File(outputDir, "unfiltered");
			else
				freqDir = new File(outputDir, (float)lowpassFreq+"hz");
			Preconditions.checkState(freqDir.exists() || freqDir.mkdir());
			BBP_Wrapper wrapper = new BBP_Wrapper(vm, Method.GP, srcFile, null, srfFile, sitesFile, freqDir);
			wrapper.setLowpassFreq(lowpassFreq);
			wrapper.setDoHF(false);
			wrapper.setDoFAS(true);
			wrapper.setDoRotD100(true);
			wrapper.setDoRotD50(false);
			wrapper.setDataOnly(true);
			
			wrapper.run();
			
			for (BBP_Site site : sites) {
				String siteName = site.getName();
				File fasFile = SpectraPlotter.findFASFile(freqDir, siteName);
				SpectraPlotter.plotFAS(fasFile, freqDir, siteName+"_fas");
				File rsRD100File = SpectraPlotter.findRotD100File(freqDir, siteName);
				SpectraPlotter.plotRotD(rsRD100File, freqDir, siteName+"_rd", true, true, null);
				
				File seisFile = SeismogramPlotter.findBBP_SeisFile(freqDir, siteName, false);
				DiscretizedFunc[] seis = SeismogramPlotter.loadBBP_Seis(seisFile);
				SeismogramPlotter.plotSeismograms(seis, "Seismograms", false, freqDir, siteName+"_seis", false, null);
			}
		}
	}

}

package scratch.kevin.bbp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CompletableFuture;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.data.function.DiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.sha.earthquake.FocalMechanism;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_Module.Method;
import scratch.kevin.bbp.BBP_Module.VelocityModel;
import scratch.kevin.bbp.BBP_SourceFile.BBP_PlanarSurface;
import scratch.kevin.simulators.ruptures.RSQSimBBP_Config;

public class BBP_LonSpanTests {

	public static void main(String[] args) throws IOException {
		File outputDir = new File("/tmp/bbp_lon_span_tests");
		Preconditions.checkState(outputDir.exists() || outputDir.mkdir());
		
		BBP_Site site1799 = new BBP_Site("site1795", new Location(0, 179.5), 760,
				RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ);
		BBP_Site siteNeg1799 = new BBP_Site("siteNeg1795", new Location(0, -179.5), 760,
				RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ);
		BBP_Site site1800 = new BBP_Site("site1800", new Location(0, 180), 760,
				RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ);
		BBP_Site site1801 = new BBP_Site("site1805", new Location(0, 180.5), 760,
				RSQSimBBP_Config.SITE_LO_PASS_FREQ, RSQSimBBP_Config.SITE_HI_PASS_FREQ);
		
		double mag = 7.5d;
		double len = 40d;
		double ddw = 10d;
		double dlen = 0.1;
		double dwid = 0.1;
		double cornerFreq = 1d;
		FocalMechanism mech = new FocalMechanism(90d, 90d, 0);
		
		BBP_SourceFile src1795 = new BBP_SourceFile(new BBP_PlanarSurface(site1799.getLoc(), len, ddw, mech), mag,
				0d, 0d, dwid, dlen, cornerFreq, 1234);
		BBP_SourceFile srcNeg1795 = new BBP_SourceFile(new BBP_PlanarSurface(siteNeg1799.getLoc(), len, ddw, mech), mag,
				0d, 0d, dwid, dlen, cornerFreq, 1234);
		BBP_SourceFile src1800 = new BBP_SourceFile(new BBP_PlanarSurface(site1800.getLoc(), len, ddw, mech), mag,
				0d, 0d, dwid, dlen, cornerFreq, 1234);
		BBP_SourceFile src1805 = new BBP_SourceFile(new BBP_PlanarSurface(site1801.getLoc(), len, ddw, mech), mag,
				0d, 0d, dwid, dlen, cornerFreq, 1234);
		
		File sitesFile = new File(outputDir, "sites.stl");
		List<BBP_Site> sites = List.of(site1799, siteNeg1799, site1800, site1801);
		BBP_Site.writeToFile(sitesFile, sites);
		List<BBP_SourceFile> sources = List.of(src1795, srcNeg1795, src1800, src1805);
		List<String> sourceNames = List.of("source1795", "sourceNeg1795", "source1800", "source1805");
		
		double[][] gm0p1s = new double[sites.size()][sources.size()];
		double[][] gm5s = new double[sites.size()][sources.size()];
		
		List<CompletableFuture<Void>> futures = new ArrayList<>();
		for (int s1=0; s1<sources.size(); s1++) {
			String srcName = sourceNames.get(s1);
			
			File srcFile = new File(outputDir, srcName+".src");
			sources.get(s1).writeToFile(srcFile);
			
			File srcOutputDir = new File(outputDir, srcName);
			Preconditions.checkState(srcOutputDir.exists() || srcOutputDir.mkdir());
			BBP_Wrapper wrapper = new BBP_Wrapper(VelocityModel.LA_BASIN_500, Method.GP, srcFile, null, null, sitesFile, srcOutputDir);
			wrapper.setDoHF(false);
			
			int sourceIndex = s1;
			futures.add(CompletableFuture.runAsync(new Runnable() {
				
				@Override
				public void run() {
					wrapper.run();
					
					for (int siteIndex=0; siteIndex<sites.size(); siteIndex++) {
						try {
							BBP_Site site = sites.get(siteIndex);
							File rd50File = SpectraPlotter.findRotD50File(srcOutputDir, site.getName());
							Preconditions.checkState(rd50File != null && rd50File.exists());
							DiscretizedFunc rd50 = SpectraPlotter.loadRotD50(rd50File);
							
							gm0p1s[sourceIndex][siteIndex] = rd50.getInterpolatedY(0.1d);
							gm5s[sourceIndex][siteIndex] = rd50.getInterpolatedY(5d);
						} catch (Exception e) {
							e.printStackTrace();
							System.exit(1);
						}
					}
				}
			}));
			try {
				Thread.sleep(2000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		System.out.println("Waiting on "+futures.size()+" futures...");
		for (CompletableFuture<Void> future : futures)
			future.join();
		System.out.println("Writing CSV");
		
		CSVFile<String> csv = new CSVFile<>(false);
		
		for (double period : new double[] {0.1d, 5d}) {
			double[][] gms = period < 1d ? gm0p1s : gm5s;
			
			if (csv.getNumRows() > 0)
				csv.addLine("");
			csv.addLine((float)period+"s RotD50 SA");
			List<String> siteHeader = new ArrayList<>();
			siteHeader.add("");
			for (BBP_Site site : sites)
				siteHeader.add(site.getName());
			csv.addLine(siteHeader);
			for (int source=0; source<sources.size(); source++) {
				List<String> line = new ArrayList<>();
				line.add(sourceNames.get(source));
				for (int site=0; site<sites.size(); site++)
					line.add((float)gms[source][site]+"");
				csv.addLine(line);
			}
		}
		
		csv.writeToFile(new File(outputDir, "gm_summary.csv"));
	}

}

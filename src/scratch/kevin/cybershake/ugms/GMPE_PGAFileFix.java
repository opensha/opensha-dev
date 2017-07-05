package scratch.kevin.cybershake.ugms;

import java.io.File;
import java.nio.ByteOrder;
import java.util.List;

import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.geo.Location;
import org.opensha.commons.util.binFile.BinaryDoubleScalarRandomAccessFile;
import org.opensha.commons.util.binFile.BinaryGeoDatasetRandomAccessFile;
import org.opensha.sha.calc.hazardMap.BinaryHazardCurveReader;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class GMPE_PGAFileFix {

	public static void main(String[] args) throws Exception {
		File mainDir = new File("/home/kevin/CyberShake/MCER/gmpe_cache_gen/");
		String dirPrefix = "2017_05_19-ucerf3_downsampled_ngaw2_binary_0.02_";
		String[] pgaNames = { "NGAWest_2014_NoIdr_MeanUCERF3_downsampled_RotD100_pga.bin",
				"NGAWest_2014_NoIdr_MeanUCERF3_downsampled_RotD100_pga_det.bin",
				"NGAWest_2014_NoIdr_MeanUCERF3_downsampled_RotD100_pga_prob.bin"};
		String mcerName = "NGAWest_2014_NoIdr_MeanUCERF3_downsampled_RotD100_mcer.bin";
		
		
		for (File dir : mainDir.listFiles()) {
			String name = dir.getName();
			if (!dir.isDirectory() || !name.startsWith(dirPrefix))
				continue;
			
			System.out.println("Processing "+name);
			
			File mcerFile = new File(dir, mcerName);
			BinaryHazardCurveReader spectrumReader = new BinaryHazardCurveReader(mcerFile.getAbsolutePath());
			List<Location> locs = Lists.newArrayList();
			ArbitrarilyDiscretizedFunc curve = spectrumReader.nextCurve();
			while (curve != null) {
				locs.add(spectrumReader.currentLocation());
				curve = spectrumReader.nextCurve();
			}
			
			for (String fileName : pgaNames) {
				System.out.println("\t"+fileName);
				File pgaFile = new File(dir, fileName);
				Preconditions.checkState(pgaFile.exists());
				File tempOut = new File(pgaFile.getAbsolutePath()+".tmp");
				double[] vals = BinaryDoubleScalarRandomAccessFile.readFile(pgaFile);
				Preconditions.checkState(vals.length == locs.size());
				BinaryGeoDatasetRandomAccessFile geoOut =
						new BinaryGeoDatasetRandomAccessFile(tempOut, ByteOrder.BIG_ENDIAN, vals.length);
				geoOut.initialize();
				for (int i=0; i<vals.length; i++)
					geoOut.write(i, locs.get(i), vals[i]);
				geoOut.close();
				Files.move(tempOut, pgaFile);
			}
		}
	}

}

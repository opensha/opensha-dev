package scratch.kevin.ucerf3.etas;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.util.FileUtils;

import scratch.UCERF3.erf.ETAS.ETAS_CatalogIO;
import scratch.UCERF3.erf.ETAS.ETAS_EqkRupture;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

public class EtasMagDepZip {
	
	private static byte[] buffer = new byte[18024];

	public static void main(String[] args) throws IOException {
		Preconditions.checkArgument(args.length == 2 || args.length == 3,
					"Usage: <results-dir> <output-file> [<min-mag>]");
		
		File resultsDir = new File(args[0]);
		Preconditions.checkArgument(resultsDir.exists() && resultsDir.isDirectory(),
				"Bad results dir: "+resultsDir.getAbsolutePath());
		
		File outputFile = new File(args[1]).getAbsoluteFile();
		Preconditions.checkArgument(outputFile.getParentFile().exists(),
				"Cannot create output: "+outputFile.getAbsolutePath());
		
		double minMag;
		if (args.length == 3)
			minMag = Double.parseDouble(args[2]);
		else
			minMag = 5;

		ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
		
		// Set the compression ratio
		out.setLevel(Deflater.DEFAULT_COMPRESSION);
		
		int count = 0;
		
		DecimalFormat df = new DecimalFormat("0.#");
		
		for (File dir : resultsDir.listFiles()) {
			if (!dir.isDirectory())
				continue;
			String name = dir.getName();
			if (!name.startsWith("sim_"))
				continue;
			
			if (count % 1000 == 0)
				System.out.println("Processing catalog "+count);
			
			String prefix = name+"/";
			
			// add the directory itself
			out.putNextEntry(new ZipEntry(prefix));
			
			addToZip(out, new File(dir, "infoString.txt"), prefix+"infoString.txt");
			
			// now sub catalog
			List<ETAS_EqkRupture> catalog;
			try {
				catalog = ETAS_CatalogIO.loadCatalog(new File(dir, "simulatedEvents.txt"), minMag);
			} catch (RuntimeException e) {
				System.out.println("Skipping "+name+": "+e.getMessage());
				continue;
			}
			File smallCat = new File(dir, "simulatedEvents_m"+df.format(minMag)+".txt");
			ETAS_CatalogIO.writeEventDataToFile(smallCat, catalog);
			
			addToZip(out, smallCat, prefix+"simulatedEvents.txt");
			
			count++;
		}
		
		// Close the ZipOutPutStream
		out.close();
	}
	
	private static void addToZip(ZipOutputStream out, File file, String name) throws IOException {
		// Add ZIP entry to output stream.
		out.putNextEntry(new ZipEntry(name));

		// Associate a file input stream for the current file
		FileInputStream in = new FileInputStream(file);

		// Transfer bytes from the current file to the ZIP file
		//out.write(buffer, 0, in.read(buffer));

		int len;
		while ((len = in.read(buffer)) > 0) {
			out.write(buffer, 0, len);
		}

		// Close the current entry
		out.closeEntry();

		// Close the current file input stream
		in.close();
	}

}

package scratch.kevin.nshm23;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.CSVFile;

import com.google.common.base.Preconditions;
import com.google.common.io.ByteStreams;
import com.google.common.io.Files;

public class SiteHazardDurationFix {
	
	public static void main(String[] args) throws ZipException, IOException {
		File input = new File(args[0]);
		Preconditions.checkState(input.exists());
		File output = new File(args[1]);
		
		double origDuration = 30d;
		double newDuration = 1d;
		
		if (input.isDirectory()) {
			Preconditions.checkState((output.exists() && output.isDirectory()) || output.mkdir());
			for (File file : input.listFiles()) {
				if (file.getName().endsWith(".csv") && (file.getName().contains("sa")
						|| file.getName().contains("pga") || file.getName().contains("pgv"))) {
					System.out.println("processing file: "+file.getName());
					CSVFile<String> csv = CSVFile.readFile(file, true);
					fixCSV(csv, origDuration, newDuration);
					// overwrite in place
					csv.writeToFile(new File(output, file.getName()));
				} else {
					System.out.println("skipping "+file.getName());
				}
			}
		} else {
			// zip file
			ZipFile zip = new ZipFile(input);
			Enumeration<? extends ZipEntry> entries = zip.entries();
			
			File tmpOutput = new File(output.getParentFile(), output.getName()+".tmp");
			ZipOutputStream zout = new ZipOutputStream(new FileOutputStream(tmpOutput));
			while (entries.hasMoreElements()) {
				ZipEntry entry = entries.nextElement();
				System.out.println("processing entry: "+entry.getName());
				if (entry.getName().endsWith(".csv") && (entry.getName().contains("sa")
						|| entry.getName().contains("pga") || entry.getName().contains("pgv"))) {
					CSVFile<String> csv = CSVFile.readStream(zip.getInputStream(entry), true);
					fixCSV(csv, origDuration, newDuration);
					
					zout.putNextEntry(new ZipEntry(entry.getName()));
					csv.writeToStream(zout);
					zout.flush();
					zout.closeEntry();
				} else {
					// copy directly
					System.out.println("copying directly");
					InputStream is = zip.getInputStream(entry);
					zout.putNextEntry(new ZipEntry(entry.getName()));
					ByteStreams.copy(new BufferedInputStream(is), zout);
					zout.flush();
					zout.closeEntry();
					is.close();
				}
			}
			zip.close();
			zout.close();
			Files.move(tmpOutput, output);
		}
	}
	
	private static void fixCSV(CSVFile<String> csv, double origDuration, double newDuration) {
		int firstDataCol= -1;
		for (int i=0; i<csv.getNumCols(); i++) {
			try {
				csv.getDouble(0, i);
				// if we got here, it worked;
				firstDataCol = i;
				break;
			} catch (NumberFormatException e) {};
		}
		Preconditions.checkState(firstDataCol > 0);
		
		for (int row=1; row<csv.getNumRows(); row++) {
			for (int col=firstDataCol; col<csv.getNumCols(); col++) {
				double origProb = csv.getDouble(row, col);
				double rate = -Math.log(1 - origProb)/origDuration;
				double newProb = 1d - Math.exp(-rate*newDuration);
				csv.set(row, col, newProb+"");
			}
		}
	}

}

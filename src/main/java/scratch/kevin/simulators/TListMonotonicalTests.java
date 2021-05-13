package scratch.kevin.simulators;

import java.io.BufferedInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.DoubleBuffer;
import java.util.Arrays;

import org.opensha.sha.simulators.parsers.RSQSimFileReader;

import com.google.common.base.Preconditions;
import com.google.common.io.LittleEndianDataInputStream;

import scratch.kevin.simulators.RSQSimCatalog.CatEnumDateComparator;
import scratch.kevin.simulators.RSQSimCatalog.Catalogs;

public class TListMonotonicalTests {

	public static void main(String[] args) throws IOException {
		Catalogs[] cats = Catalogs.values();
		Arrays.sort(cats, new CatEnumDateComparator());
		
		File baseDir = new File("/data/kevin/simulators/catalogs");
		
		for (Catalogs cat : cats) {
			System.out.println("Processing "+cat.name());
			
			File catDir = cat.instance(baseDir).getCatalogDir();
			File tListFile;
			try {
				tListFile = RSQSimFileReader.findByExt(catDir, "tList");
			} catch (Exception e) {
				System.out.println("Skipping");
				continue;
			}
			
			LittleEndianDataInputStream din = new LittleEndianDataInputStream(new BufferedInputStream(new FileInputStream(tListFile)));
			
			int numViolations = 0;
			
			double prevTime = Double.NaN;
			
			int index = 0;
			while (true) {
				double time;
				try {
					time = din.readDouble();
				} catch (EOFException e) {
					break;
				}
				Preconditions.checkState(Double.isFinite(time), "Bad time: %s", time);
				
				if (!Double.isNaN(prevTime)) {
					if ((float)time < (float)prevTime) {
						if (numViolations < 3) {
							System.out.println("Bad time detected!");
							System.out.println("\tIndex "+(index - 1)+": "+(float)prevTime);
							System.out.println("\tIndex "+index+": "+(float)time);
						} else if (numViolations == 3) {
							System.out.println("Suppressing future violations");
						}
						numViolations++;
					}
				}
				
				prevTime = time;
				index++;
			}
			
			double percent = 100d*(double)numViolations/(double)(index-1);
			
			din.close();
			
			System.out.println(cat.name()+" had "+numViolations+" violations ("+(float)percent+" %)");
		}
	}

}

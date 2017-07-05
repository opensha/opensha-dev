package scratch.kevin;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;

import org.apache.commons.io.FileUtils;

import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;

public class UCVM_BinFileMaker {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// this generates binary files from the UCVM basin depth text files
		
		File dir = new File("/home/kevin/OpenSHA/ucvm/ucvm_basin_maps_20120521/basins");
		
		byte[] recordBuffer = new byte[4];
		ByteBuffer record = ByteBuffer.wrap(recordBuffer);
		record.order(ByteOrder.LITTLE_ENDIAN);
		
		FloatBuffer fbuff = record.asFloatBuffer();
		
		for (File file : dir.listFiles()) {
			if (file.isDirectory())
				continue;
			String name = file.getName();
			if (!name.toLowerCase().endsWith(".txt"))
				continue;
			
			name = name.substring(0, name.toLowerCase().indexOf(".txt"));
			
			OutputStream firstFW = new FileOutputStream(new File(dir, name+"_first.bin"));
			firstFW = new BufferedOutputStream(firstFW);
			OutputStream lastFW = new FileOutputStream(new File(dir, name+"_last.bin"));
			lastFW = new BufferedOutputStream(lastFW);
			
			for (String line : FileUtils.readLines(file)) {
				line = line.trim();
				if (line.isEmpty())
					continue;
				String[] strs =
						Iterables.toArray(Splitter.on(' ').omitEmptyStrings().split(line),String.class);
				
				float first = Float.parseFloat(strs[2]);
				float last = Float.parseFloat(strs[3]);
				
				fbuff.put(0, first);
				firstFW.write(recordBuffer);
				
				fbuff.put(0, last);
				lastFW.write(recordBuffer);
			}
			
			firstFW.close();
			lastFW.close();
		}
	}

}

package scratch.kevin.ucerf3.eal;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;

import com.google.common.base.Joiner;

public class GridFileDebug {

	public static void main(String[] args) throws IOException {
		DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(
				new File("/tmp/cfss_branches/FM3_1_ZENGBB_HB08_DsrTap_CharConst_M5Rate6.5"
						+ "_MMaxOff7.6_NoFix_SpatSeisU3_grid_source.bin"))));
		
		FileWriter fw = new FileWriter(new File("/tmp/grid.txt"));
		int numFuncs = in.readInt();
		for (int i=0; i<numFuncs; i++) {
			int arraySize = in.readInt();

			String str;
			if (i % 2 == 0) {
				fw.write("\nFunc "+i/2+"\n");
				str = "X:";
			} else {
				str = "Y:";
			}
			for (int j=0; j<arraySize; j++)
				str += " "+(float)in.readDouble();
			
			fw.write("\t"+str+"\n");
		}
		fw.close();
	}

}

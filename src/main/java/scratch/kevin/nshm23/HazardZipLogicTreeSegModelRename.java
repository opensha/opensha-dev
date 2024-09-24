package scratch.kevin.nshm23;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

import org.opensha.commons.util.modules.ModuleArchiveInput;
import org.opensha.commons.util.modules.ModuleArchiveOutput;
import org.opensha.sha.earthquake.faultSysSolution.modules.AbstractLogicTreeModule;

public class HazardZipLogicTreeSegModelRename {
	
	public static void main(String[] args) throws IOException {
		File inputFile = new File(args[0]);
		File outputFile = new File(args[1]);
		
		ModuleArchiveInput.ApacheZipFileInput input = new ModuleArchiveInput.ApacheZipFileInput(inputFile);
		ModuleArchiveOutput.ApacheZipFileOutput output = new ModuleArchiveOutput.ApacheZipFileOutput(outputFile);
		
		for (String entry : input.getEntries()) {
			if (entry.equals(AbstractLogicTreeModule.LOGIC_TREE_FILE_NAME) || entry.equals(AbstractLogicTreeModule.LOGIC_TREE_MAPPINGS_FILE_NAME)) {
				System.out.println("Translating "+entry);
				InputStream is = input.getInputStream(entry);
				BufferedReader bRead = new BufferedReader(new InputStreamReader(is));
				
				output.putNextEntry(entry);
				
				BufferedWriter bWrite = new BufferedWriter(new OutputStreamWriter(output.getOutputStream()));
				
				String line;
				while ((line = bRead.readLine()) != null) {
					bWrite.write(replace(line));
					bWrite.write('\n');
				}
				
				bWrite.flush();
				output.closeEntry();
				bRead.close();
			} else {
				String modName = replace(entry);
//				System.out.println("Copying "+entry+" to "+modName);
				output.transferFrom(input, entry, modName);
			}
		}
		
		output.close();
		input.close();
	}
	
	private static String replace(String str) {
		return str.replaceAll("HighSeg", "High").replaceAll("MidSeg", "Middle").replaceAll("LowSeg", "Low");
	}

}

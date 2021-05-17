package scratch.kevin;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.opensha.commons.util.FileUtils;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.io.Files;

public class ShootoutRandomizer {

	public static void main(String[] args) throws IOException {
//		List<String> pres = Lists.newArrayList("312", "focusrite", "great_river", "soundcraft");
//		List<String> prefixes = Lists.newArrayList("bounce_acoustic", "bounce_vocal");
		
		List<String> pres = Lists.newArrayList("focusrite", "rosetta", "mini_me");
		List<String> prefixes = Lists.newArrayList("mix", "guitar", "vox");
		
//		List<String> participants = Lists.newArrayList("kevin", "nate", "sean", "kris", "mark", "test");
		List<String> participants = Lists.newArrayList("kevin", "nate", "extra1", "extra2");
		
		List<String> randNames = Lists.newArrayList("A", "B", "C", "D");
		
//		File dir = new File("/tmp/312 Build Shootout");
		File dir = new File("/home/kevin/Documents/studio/converter_shootout");
		
		for (String participant : participants) {
			File pDir = new File(dir, participant);
			if (!pDir.exists())
				pDir.mkdir();
			
			List<String> pFileNames = Lists.newArrayList();
			
			FileWriter fw = new FileWriter(new File(dir, participant+"_key.txt"));
			fw.write("Key for "+participant+"\n");
			
			for (String prefix : prefixes) {
				fw.write("\n"+prefix+"\n");
				List<String> myPres = Lists.newArrayList(pres);
				Collections.shuffle(myPres);
				for (int i = 0; i < myPres.size(); i++) {
					String r = randNames.get(i);
					String pre = myPres.get(i);
//					File wavFile = new File(dir, prefix+"_"+pre+".wav");
					File wavFile = new File(dir, pre+"_"+prefix+".wav");
					Preconditions.checkState(wavFile.exists(), "file not found: "+wavFile.getAbsolutePath());
					File outFile = new File(pDir, prefix+"_"+r+".wav");
					Files.copy(wavFile, outFile);
					pFileNames.add(outFile.getName());
					fw.write("\t"+r+": "+pre+"\n");
				}
			}
			fw.close();
			FileUtils.createZipFile(new File(dir, participant+".zip").getAbsolutePath(), pDir.getAbsolutePath(), pFileNames);
		}
	}

}

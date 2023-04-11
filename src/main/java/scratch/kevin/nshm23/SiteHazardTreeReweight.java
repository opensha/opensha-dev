package scratch.kevin.nshm23;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.opensha.commons.data.CSVFile;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeNode;
import org.opensha.sha.earthquake.faultSysSolution.modules.AbstractLogicTreeModule;

import com.google.common.base.Preconditions;
import com.google.common.io.ByteStreams;
import com.google.common.io.Files;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class SiteHazardTreeReweight {

	public static void main(String[] args) throws ZipException, IOException {
		File input = new File(args[0]);
		Preconditions.checkState(input.exists());
		File output = new File(args[1]);
		File ltFile = new File(args[2]);
		
		LogicTree<?> tree = LogicTree.read(ltFile);
		
		Map<String, LogicTreeBranch<?>> branchStringsMap = new HashMap<>();
		
		if (input.isDirectory()) {
			Preconditions.checkState((output.exists() && output.isDirectory()) || output.mkdir());
			for (File file : input.listFiles()) {
				if (file.getName().endsWith(".csv") && (file.getName().contains("sa")
						|| file.getName().contains("pga") || file.getName().contains("pgv"))) {
					System.out.println("processing file: "+file.getName());
					CSVFile<String> csv = CSVFile.readFile(file, true);
					fixCSV(csv, tree, branchStringsMap);
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
					fixCSV(csv, tree, branchStringsMap);
					
					zout.putNextEntry(new ZipEntry(entry.getName()));
					csv.writeToStream(zout);
					zout.flush();
					zout.closeEntry();
				} else if (entry.getName().equals(AbstractLogicTreeModule.LOGIC_TREE_FILE_NAME)) {
					entry = new ZipEntry(AbstractLogicTreeModule.LOGIC_TREE_FILE_NAME);
					zout.putNextEntry(entry);
					Gson gson = new GsonBuilder().setPrettyPrinting()
							.registerTypeAdapter(LogicTree.class, new LogicTree.Adapter<>()).create();
					BufferedOutputStream out = new BufferedOutputStream(zout);
					BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(out));
					gson.toJson(tree, LogicTree.class, writer);
					writer.flush();
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
	
	private static String BR_KEY_SEP = "||";
	
	private static void fixCSV(CSVFile<String> csv, LogicTree<?> tree, Map<String, LogicTreeBranch<?>> branchStringsMap) {
		int branchIndCol = 1;
		int branchWeightCol = 2;
		int firstBranchCol = 3;
		Preconditions.checkState(tree.size() == csv.getNumRows()-1);
		
		for (int row=1; row<csv.getNumRows(); row++) {
			int index = csv.getInt(row, branchIndCol);
			LogicTreeBranch<?> branch = tree.getBranch(index);
			boolean match = true;
			String branchKey = null;
			for (int i=0; i<branch.size(); i++) {
				String name = branch.getValue(i).getShortName();
				String csvName = csv.get(row, firstBranchCol+i);
				if (branchKey == null)
					branchKey = "";
				else
					branchKey += BR_KEY_SEP;
				branchKey += csvName;
				if (!name.equals(csvName)) {
					match = false;
				}
			}
			if (!match) {
				// need to manually match
				if (branchStringsMap.isEmpty()) {
					// need to init
					for (LogicTreeBranch<?> oBranch : tree) {
						String oBranchKey = null;
						for (LogicTreeNode node : oBranch) {
							if (oBranchKey == null)
								oBranchKey = "";
							else
								oBranchKey += BR_KEY_SEP;
							oBranchKey += node.getShortName();
						}
						Preconditions.checkState(!branchStringsMap.containsKey(oBranchKey));
						branchStringsMap.put(oBranchKey, oBranch);
					}
				}
				branch = branchStringsMap.get(branchKey);
				Preconditions.checkNotNull(branch);
			}
			csv.set(row, branchWeightCol, tree.getBranchWeight(branch)+"");
		}
	}

}

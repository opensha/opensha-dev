package scratch.kevin.nshm23;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.commons.compress.archivers.zip.ZipArchiveEntry;
import org.apache.commons.compress.archivers.zip.ZipArchiveEntryPredicate;
import org.apache.commons.compress.archivers.zip.ZipArchiveOutputStream;
import org.apache.commons.compress.archivers.zip.ZipFile;
import org.opensha.commons.logicTree.BranchWeightProvider;
import org.opensha.commons.logicTree.LogicTree;
import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.util.ExceptionUtils;

import com.google.common.io.Files;

public class UpdateLogicTreeWeights {

	public static void main(String[] args) throws IOException {
		if (args.length == 0 || args.length > 2) {
			System.err.println("USAGE: <logic-tree.json OR results.zip> [<output-file>]");
			System.exit(0);
		}
		
		File file = new File(args[0]);
		File outputFile = args.length == 1 ? file : new File(args[1]);
		
		BranchWeightProvider prov = new BranchWeightProvider.CurrentWeights();
		
		if (file.getName().endsWith(".zip")) {
			// zip file
			
			File tmpOut = new File(outputFile.getAbsolutePath()+".tmp");
			ZipArchiveOutputStream zout = new ZipArchiveOutputStream(tmpOut);
			
			ZipFile zip = new ZipFile(file);
			
			zip.copyRawEntries(zout, new ZipArchiveEntryPredicate() {
				
				@Override
				public boolean test(ZipArchiveEntry zipArchiveEntry) {
					if (zipArchiveEntry.getName().endsWith("logic_tree.json")) {
						// update it
						try {
							System.out.println("Updaing logic tree archive: "+zipArchiveEntry.getName());
							BufferedReader reader = new BufferedReader(new InputStreamReader(zip.getInputStream(zipArchiveEntry)));
							LogicTree<?> tree = LogicTree.read(reader);
							updateOrigWeights(tree, prov);
							
							zout.putArchiveEntry(new ZipArchiveEntry(zipArchiveEntry.getName()));
							BufferedOutputStream bStream = new BufferedOutputStream(zout);
							tree.writeToStream(bStream);
							bStream.flush();
							zout.flush();
							zout.closeArchiveEntry();
						} catch (IOException e) {
							throw ExceptionUtils.asRuntimeException(e);
						}
						return false;
					}
					return true;
				}
			});
			
			zip.close();
			zout.close();
			Files.move(tmpOut, outputFile);
		} else {
			// assume JSON
			LogicTree<?> tree = LogicTree.read(file);
			updateOrigWeights(tree, prov);
			tree.write(outputFile);
		}
	}
	
	private static void updateOrigWeights(LogicTree<?> tree, BranchWeightProvider prov) {
		for (LogicTreeBranch<?> branch : tree)
			branch.setOrigBranchWeight(prov.getWeight(branch));
	}

}

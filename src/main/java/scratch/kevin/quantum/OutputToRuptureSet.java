package scratch.kevin.quantum;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.opensha.sha.earthquake.faultSysSolution.FaultSystemRupSet;

import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.stream.JsonReader;

public class OutputToRuptureSet {

	public static void main(String[] args) throws IOException {
		File origDir = new File("/home/kevin/markdown/inversions/2024_11_06-quantum_test_rup_set_problem-213_sects");
		File origExhaustiveRupSetFile = new File(origDir, "orig_rup_set.zip");
		File origThinnedRupSetFile = new File(origDir, "orig_rup_set.zip");
		
		File quantumOutputFile = new File("/home/kevin/OpenSHA/quantum/rapture_gen/embedding213.json");
		
		FaultSystemRupSet origThinnedRupSet = FaultSystemRupSet.load(origThinnedRupSetFile);
		FaultSystemRupSet origExhaustiveRupSet = FaultSystemRupSet.load(origExhaustiveRupSetFile);
		
		Gson gson = new GsonBuilder().create();
		
		JsonReader in = gson.newJsonReader(new BufferedReader(new FileReader(quantumOutputFile)));
		
		int numSects = origExhaustiveRupSet.getNumSections();
		
		int[] sectCounts = new int[numSects];
		
		List<List<Integer>> rupSectsList = new ArrayList<>();
		
		in.beginObject();
		
		while (in.hasNext()) {
			System.out.println("Reading "+in.nextName());
			in.beginArray();
			List<Integer> sectIDs = new ArrayList<>();
			while (in.hasNext()) {
				int id = in.nextInt();
				Preconditions.checkState(id >= 0 && id < numSects, "Bad id=%s, numSects=%s", id, numSects);
				Preconditions.checkState(!sectIDs.contains(id), "Rupture contains a duplicate ID: %s", id);
				sectCounts[id]++;
				sectIDs.add(id);
			}
			in.endArray();
			rupSectsList.add(sectIDs);
		}
		
		in.endObject();
		
		in.close();
	}

}

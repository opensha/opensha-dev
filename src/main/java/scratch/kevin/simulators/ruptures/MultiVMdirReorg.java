package scratch.kevin.simulators.ruptures;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;

import org.opensha.commons.util.FileNameComparator;

import com.google.common.base.Preconditions;

import scratch.kevin.bbp.BBP_Module.VelocityModel;

public class MultiVMdirReorg {

	public static void main(String[] args) throws IOException {
		String[] wildcards = {
				"gmpe_bbp_comparisons_",
				"event_",
				"gmpe_bbp_rg_comparisons_",
				"catalog_rotd_ratio_comparisons",
				"source_site_comparisons",
				"bbp_part_b",
				"rotated_ruptures_"
		};
		
		File gitDir = new File("/home/kevin/markdown/rsqsim-analysis/catalogs");
		
		VelocityModel defaultVM = VelocityModel.LA_BASIN_863;
		
		boolean dryRun = false;
		
		File[] dirs = gitDir.listFiles();
		Arrays.sort(dirs, new FileNameComparator());
		
		for (File dir : dirs) {
			if (!dir.isDirectory())
				continue;
			System.out.println("Processing "+dir.getName());
			
			for (VelocityModel vm : VelocityModel.values()) {
				File oldDir = new File(dir, "bbp_vm_"+vm.name());
				if (oldDir.exists()) {
					File dest = new File(dir, "bbp_"+vm.name());
					if (!dryRun)
						Files.move(oldDir.toPath(), dest.toPath());
				}
			}
			
			for (File subDir : dir.listFiles()) {
				if (!subDir.isDirectory())
					continue;
				
				String name = subDir.getName();
				
				boolean match = false;
				for (String wildcard : wildcards)
					if (name.contains(wildcard))
						match = true;
				
				if (match) {
					VelocityModel myVM = null;
					for (VelocityModel vm : VelocityModel.values()) {
						if (name.contains(vm.name()))
							myVM = vm;
					}
					if (myVM == null)
						myVM = defaultVM;
					File bbpDir = new File(dir, "bbp_"+myVM.name());
					Preconditions.checkState(dryRun || bbpDir.exists() || bbpDir.mkdir());
					String destName = name.replaceAll("_vm"+myVM.name(), "");
					File destDir = new File(bbpDir, destName);
					System.out.println("\t"+name+"\t=>\t"+destDir.getAbsolutePath().replaceAll(dir.getAbsolutePath(), ""));
					if (!dryRun)
						Files.move(subDir.toPath(), destDir.toPath());
				} else {
					System.out.println("\t"+name+"\t\t(do nothing)");
				}
			}
		}
	}

}

package scratch.nshm23.logicTree;

import java.util.ArrayList;
import java.util.List;

import org.opensha.commons.logicTree.LogicTreeBranch;
import org.opensha.commons.logicTree.LogicTreeLevel;
import org.opensha.commons.logicTree.LogicTreeNode;

import scratch.UCERF3.enumTreeBranches.DeformationModels;
import scratch.UCERF3.enumTreeBranches.FaultModels;
import scratch.UCERF3.enumTreeBranches.ScalingRelationships;
import scratch.UCERF3.enumTreeBranches.SlipAlongRuptureModels;

public class DraftNSHM23LogicTreeBranch extends LogicTreeBranch<LogicTreeNode> {
	
	public static List<LogicTreeLevel<? extends LogicTreeNode>> levels;
	
	static {
		levels = new ArrayList<>();
		levels.add(LogicTreeLevel.forEnum(FaultModels.class, "Fault Model", "FM"));
		levels.add(LogicTreeLevel.forEnum(DeformationModels.class, "Deformation Model", "DM"));
		levels.add(LogicTreeLevel.forEnum(ScalingRelationships.class, "Scaling Relationship", "Scale"));
		levels.add(LogicTreeLevel.forEnum(SlipAlongRuptureModels.class, "Slip Along Rupture", "SlipAlong"));
		levels.add(LogicTreeLevel.forEnum(SupraSeisBValue.class, "Supra-Seis b-value", "SupraB"));
		levels.add(LogicTreeLevel.forEnum(SubSectConstraintModel.class, "Sub-Sect Constraint Model", "SectConstr"));
	}
	
	private DraftNSHM23LogicTreeBranch() {
		super(levels);
	}
	
	public DraftNSHM23LogicTreeBranch(List<LogicTreeNode> values) {
		super(levels, values);
	}

}

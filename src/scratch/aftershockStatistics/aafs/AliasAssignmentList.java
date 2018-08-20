package scratch.aftershockStatistics.aafs;

import java.util.Collection;
import java.util.Set;
import java.util.HashSet;
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Queue;
import java.util.Deque;
import java.util.ArrayDeque;
import java.util.Comparator;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;


/**
 * A list of timelines, containing the assignment of Comcat IDs to each timeline.
 * Author: Michael Barall 07/01/2018.
 *
 * Given any pair of timelines in the list, the set of current Comcat IDs for the
 * two timelines is required to be disjoint.  There is no such requirement for removed IDs.
 */
public class AliasAssignmentList {


	//----- Values -----

	// The list of assignments.

	private ArrayList<AliasAssignment> assignments;

	// The map keys are all the timeline IDs that are currently assigned to any timeline.
	// The map value is the timeline that contains the ID as its timeline ID.
	// Note that we require only one timeline can have a given timeline ID.
	// Note this map may not refer to all assignments, since it is possible for an
	// assignment to have null timeline ID.

	private HashMap<String, AliasAssignment> all_timeline_ids;

	// The map keys are all the Comcat IDs that are currently assigned to any timeline.
	// The map value is the timeline that contains the ID as a current ID.
	// Note that we require only one timeline can contain a Comcat ID as a current ID.

	private HashMap<String, AliasAssignment> all_comcat_ids;

	// The Comcat IDs that were previously assigned to any timeline, but have been removed,
	// and are not currently assigned to any timeline.
	// Note that all_removed_ids and all_comcat_ids.keySet() are disjoint sets.
	// Their union is the set of all Comcat IDs that are mentioned in any timeline.

	private HashSet<String> all_removed_ids;

	// A time stamp, in milliseconds since the epoch, or 0L if none.

	long family_time;

	// The Comcat IDs that are known to be absent from the underlying data store.

	private HashSet<String> all_absent_ids;


	//----- Construction -----

	// The default constructor creates an empty assignment list.

	public AliasAssignmentList () {
		this.assignments = new ArrayList<AliasAssignment>();
		this.all_timeline_ids = new HashMap<String, AliasAssignment>();
		this.all_comcat_ids = new HashMap<String, AliasAssignment>();
		this.all_removed_ids = new HashSet<String>();
		this.family_time = 0L;
		this.all_absent_ids = new HashSet<String>();
	}

	// Clear the assignment list.

	public void clear () {
		this.assignments = new ArrayList<AliasAssignment>();
		this.all_timeline_ids = new HashMap<String, AliasAssignment>();
		this.all_comcat_ids = new HashMap<String, AliasAssignment>();
		this.all_removed_ids = new HashSet<String>();
		this.family_time = 0L;
		this.all_absent_ids = new HashSet<String>();
	}


	//----- Access -----

	// Get the number of assignments.

	public int get_assignment_count () {
		return assignments.size();
	}


	// Get the i-th assignment.

	public AliasAssignment get_assignment (int i) {
		return assignments.get(i);
	}


	// Get an iterator over assignments.

	public Iterator<AliasAssignment> get_assigment_iterator () {
		return assignments.iterator();
	}


	// Get the time stamp.

	public long get_family_time () {
		return family_time;
	}


	// Set the time stamp.

	public void set_family_time (long the_family_time) {
		family_time = the_family_time;
		return;
	}


	// Get the assignment for a given timeline ID.
	// Returns null if none found.

	public AliasAssignment get_assignment_for_timeline_id (String timeline_id) {
		return all_timeline_ids.get (timeline_id);
	}


	// Get the assignment for a given Comcat ID.
	// Returns null if none found.
	// Note: This considers only currently assigned Comcat IDs, not previously
	// assigned Comcat IDs, which implies there can be only one assignment.

	public AliasAssignment get_assignment_for_comcat_id (String comcat_id) {
		return all_comcat_ids.get (comcat_id);
	}


	// Get the assignment for a given primary ID.
	// Returns null if none found.

	public AliasAssignment get_assignment_for_primary_id (String primary_id) {
		AliasAssignment assignment = all_comcat_ids.get (primary_id);
		if (assignment != null) {
			if (!( primary_id.equals (assignment.get_primary_id()) )) {		// equals is always false if its argument is null
				assignment = null;
			}
		}
		return assignment;
	}


//	// Add an assignment to the list.
//	// An exception is thrown if any current Comcat ID is already a current Comcat ID for another assignment.
//	// Note that the assignment need not have a timeline ID defined.
//
//	public void add_assignment (AliasAssignment assignment) {
//	
//		// Add it to the List
//
//		assignments.add (assignment);
//
//		// Add its IDs to the tables
//
//		assignment.accumulate_all_ids (all_comcat_ids, all_removed_ids);
//
//		// If it has a timeline ID, add it to the tables
//
//		String timeline_id = assignment.get_timeline_id();
//		if (timeline_id != null) {
//			if (all_timeline_ids.put (timeline_id, assignment) != null) {
//				throw new IllegalArgumentException("AliasAssignmentList.add_assignment: Duplicate timeline ID supplied : " + timeline_id);
//			}
//		}
//
//		return;
//	}


	// Add an assignment to the list.
	// Note that the assignment need not have a timeline ID defined.
	// An exception is thrown if a duplicate Comcat ID or timeline ID is encountered.
	// Note: It is not an error for the new assignment to have a Comcat ID that is already in
	// all_removed_ids, or to have a removed ID that is already in all_comcat_ids or all_removed_ids.

	public void add_assignment (AliasAssignment assignment) {
	
		// Add it to the List

		assignments.add (assignment);

		// Add Comcat IDs to all_comcat_ids, and remove them from all_removed_ids

		for (Iterator<String> it = assignment.get_comcat_id_iterator(); it.hasNext(); ) {
			String comcat_id = it.next();
			if (all_comcat_ids.put (comcat_id, assignment) != null) {
				throw new IllegalArgumentException("AliasAssignmentList.add_assignment: Duplicate Comcat ID supplied : " + comcat_id);
			}
			all_removed_ids.remove (comcat_id);
		}

		// Add removed IDs to all_removed_ids, if they are not already Comcat IDs

		for (Iterator<String> it = assignment.get_removed_id_iterator(); it.hasNext(); ) {
			String removed_id = it.next();
			if (!( all_comcat_ids.containsKey (removed_id) )) {
				all_removed_ids.add (removed_id);
			}
		}

		// If it has a timeline ID, add it to the tables

		String timeline_id = assignment.get_timeline_id();
		if (timeline_id != null) {
			if (all_timeline_ids.put (timeline_id, assignment) != null) {
				throw new IllegalArgumentException("AliasAssignmentList.add_assignment: Duplicate timeline ID supplied : " + timeline_id);
			}
		}

		return;
	}


	// Get the total number of Comcat IDs.

	public int get_all_comcat_ids_count () {
		return all_comcat_ids.size();
	}


	// Get an iterator over Comcat IDs.

	public Iterator<String> get_all_comcat_id_iterator () {
		return all_comcat_ids.keySet().iterator();
	}


	// Return true if the given ID is in the list of Comcat IDs.

	public boolean contains_comcat_id (String id) {
		return all_comcat_ids.containsKey(id);
	}


	// Get the total number of removed IDs.

	public int get_all_removed_ids_count () {
		return all_removed_ids.size();
	}


	// Get an iterator over removed IDs.

	public Iterator<String> get_all_removed_id_iterator () {
		return all_removed_ids.iterator();
	}


	// Return true if the given ID is in the list of removed IDs.

	public boolean contains_removed_id (String id) {
		return all_removed_ids.contains(id);
	}


	// Get the total number of removed IDs, summed over all assignments.
	// Note this can be larger than get_all_removed_ids_count() because removed IDs
	// can be duplicated among assignments.

	public int get_total_removed_ids_count () {
		int result = 0;
		for (AliasAssignment assignment : assignments) {
			result += assignment.get_removed_id_count();
		}
		return result;
	}


	// Get all the removed IDs into an array range.
	// The range extends from index lo (inclusive) to index hi (exclusive),
	// where hi = lo + get_all_removed_ids_count().
	// If lo == hi then the range is empty, in which case nothing is retrieved.
	// The array must be non-null, and must satisfy 0 <= lo <= hi <= idlist.length.
	// The return value is hi.

	public int get_all_removed_ids_as_array (String[] idlist, int lo) {

		// Check arguments

		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignmentList.get_all_removed_ids_as_array: No array supplied");
		}

		int hi = lo + all_removed_ids.size();

		if (!( 0 <= lo && lo <= hi && hi <= idlist.length )) {
			throw new IllegalArgumentException("AliasAssignmentList.get_all_removed_ids_as_array: Invalid range : lo = " + lo + ", hi = " + hi + ", length = " + idlist.length);
		}

		// Insert the removed IDs

		int i = lo;
		for (String id : all_removed_ids) {
			idlist[i++] = id;
		}

		return hi;
	}


	// Get all the removed IDs into an array.
	// The length of the returned array equals get_all_removed_ids_count().

	public String[] get_all_removed_ids_as_array () {
		String[] idlist = new String[all_removed_ids.size()];
		get_all_removed_ids_as_array (idlist, 0);
		return idlist;
	}


	// Return true if the given id is known (belongs to either all_comcat_ids or all_removed_ids).

	public boolean is_known_id (String id) {
		return all_comcat_ids.containsKey (id) || all_removed_ids.contains (id);
	}


	// Add all the assignments in the other list to this list.
	// Also set the time stamp to the maximum of this and the other time stamp.
	// An exception is thrown if a duplicate Comcat ID or timeline ID is encountered.
	// Note: It is not an error for a new assignment to have a Comcat ID that is already in
	// all_removed_ids, or to have a removed ID that is already in all_comcat_ids or all_removed_ids.

	public void merge_from (AliasAssignmentList other) {

		// Loop over assignments in the other list, and add to this list

		for (AliasAssignment assignment : other.assignments) {
			add_assignment (assignment);
		}

		// Keep the maximum time stamp

		if (family_time < other.family_time) {
			family_time = other.family_time;
		}

		return;
	}


	// toString - Convert to string.

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();
		result.append ("AliasAssignmentList:" + "\n");
		result.append ("\tfamily_time: " + family_time + "\n");

		int i = 0;
		for (AliasAssignment assignment : assignments) {
			result.append (assignment.toString(i) + "\n");
			++i;
		}

		result.append ("\tall_absent_ids: " + Arrays.toString (get_all_absent_ids_as_array()));

		return result.toString();
	}


	//----- Absent IDs -----


	// Get the total number of absent IDs.

	public int get_all_absent_ids_count () {
		return all_absent_ids.size();
	}


	// Get an iterator over absent IDs.

	public Iterator<String> get_all_absent_id_iterator () {
		return all_absent_ids.iterator();
	}


	// Return true if the given ID is in the list of absent IDs.

	public boolean contains_absent_id (String id) {
		return all_absent_ids.contains(id);
	}


	// Get all the absent IDs into an array range.
	// The range extends from index lo (inclusive) to index hi (exclusive),
	// where hi = lo + get_all_absent_ids_count().
	// If lo == hi then the range is empty, in which case nothing is retrieved.
	// The array must be non-null, and must satisfy 0 <= lo <= hi <= idlist.length.
	// The return value is hi.

	public int get_all_absent_ids_as_array (String[] idlist, int lo) {

		// Check arguments

		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignmentList.get_all_absent_ids_as_array: No array supplied");
		}

		int hi = lo + all_absent_ids.size();

		if (!( 0 <= lo && lo <= hi && hi <= idlist.length )) {
			throw new IllegalArgumentException("AliasAssignmentList.get_all_absent_ids_as_array: Invalid range : lo = " + lo + ", hi = " + hi + ", length = " + idlist.length);
		}

		// Insert the absent IDs

		int i = lo;
		for (String id : all_absent_ids) {
			idlist[i++] = id;
		}

		return hi;
	}


	// Get all the absent IDs into an array.
	// The length of the returned array equals get_all_absent_ids_count().

	public String[] get_all_absent_ids_as_array () {
		String[] idlist = new String[all_absent_ids.size()];
		get_all_absent_ids_as_array (idlist, 0);
		return idlist;
	}


	// Return true if any ID in the assignment is in the list of absent IDs.

	public boolean contains_any_absent_id (AliasAssignment assignment) {

		// Check Comcat IDs

		for (Iterator<String> it = assignment.get_comcat_id_iterator(); it.hasNext(); ) {
			String comcat_id = it.next();
			if (all_absent_ids.contains (comcat_id)) {
				return true;
			}
		}

		// Add removed IDs to all_removed_ids, if they are not already Comcat IDs

		for (Iterator<String> it = assignment.get_removed_id_iterator(); it.hasNext(); ) {
			String removed_id = it.next();
			if (all_absent_ids.contains (removed_id)) {
				return true;
			}
		}

		return false;
	}


	// Return true if the given id is known (belongs to either all_comcat_ids or all_removed_ids) or absent.

	public boolean is_known_or_absent_id (String id) {
		return all_comcat_ids.containsKey (id) || all_removed_ids.contains (id) || all_absent_ids.contains (id);
	}


	// Add the given ID to the set of absent IDs.
	// Note: No checking is done, and duplicate IDs are ignored.

	public void add_absent_id (String id) {
		all_absent_ids.add (id);
		return;
	}


	//----- Matching -----


	// Sort the assignments according to timeline ID.
	// An exception is thrown if any assignment does not have a timeline ID.
	// Note: Assignment lists should be sorted before starting to match, to ensure repeatable results.

	public void sort_by_timeline_id () {

		// A comparator that compares based on timeline ID

		Comparator<AliasAssignment> comp = new Comparator<AliasAssignment>() {
			@Override public int compare (AliasAssignment o1, AliasAssignment o2) {
				return o1.compare_timeline_id_to (o2);
			}
		};

		// Sort the list if it has more than one element

		if (assignments.size() > 1) {
			assignments.sort (comp);
		}
	
		return;
	}


	// Sort the assignments according to primary ID.
	// An exception is thrown if any assignment does not have a primary ID.
	// Note: Assignment lists should be sorted before starting to match, to ensure repeatable results.

	public void sort_by_primary_id () {

		// A comparator that compares based on primary ID

		Comparator<AliasAssignment> comp = new Comparator<AliasAssignment>() {
			@Override public int compare (AliasAssignment o1, AliasAssignment o2) {
				return o1.compare_primary_id_to (o2);
			}
		};

		// Sort the list if it has more than one element

		if (assignments.size() > 1) {
			assignments.sort (comp);
		}
	
		return;
	}


	// Match a successor to a predecessor.
	// Parameters:
	//  predecessor = An assignment that is NOT in this list.  It must have a timeline ID.
	//                It must not currently have a successor.
	//  successor = An assignment in this list.  It must not currently have a predecessor.
	// The predecessor's timeline ID is assigned to the successor.  In addition, all of the
	// predecessor's Comcat IDs and removed IDs are added to the successor's removed IDs,
	// excluding those IDs that are already known to the successor.

	public void match_successor (AliasAssignment predecessor, AliasAssignment successor) {

		// Check that the successor does not already have a predecessor

		if (successor.has_match()) {
			throw new IllegalArgumentException ("AliasAssignmentList.match_successor: Successor already has a predecessor");
		}

		// Check that the predecessor does not already have a successor

		if (predecessor.has_match()) {
			throw new IllegalArgumentException ("AliasAssignmentList.match_successor: Predecessor already has a successor");
		}

		// Link the predecessor and successor

		predecessor.set_match_assignment (successor);
		successor.set_match_assignment (predecessor);

		// Get the predecessor's timeline ID and assign it to the successor

		String timeline_id = predecessor.get_timeline_id();
		if (!( timeline_id != null )) {
			throw new IllegalArgumentException ("AliasAssignmentList.match_successor: Predecessor has a null timeline ID");
		}

		successor.set_timeline_id (timeline_id);

		if (all_timeline_ids.put (timeline_id, successor) != null) {
			throw new IllegalArgumentException ("AliasAssignmentList.match_successor: Predecessor has a duplicate timeline ID: " + timeline_id);
		}

		// Take all the predecessor's Comcat IDs and removed IDs, and make them removed IDs in
		// the successor if the successor does not already know them

		successor.add_removed_ids_from_predecessor (predecessor);

		// Add the successor's removed IDs to all_removed_ids, if they are not already Comcat IDs

		for (Iterator<String> it = successor.get_removed_id_iterator(); it.hasNext(); ) {
			String removed_id = it.next();
			if (!( all_comcat_ids.containsKey (removed_id) )) {
				all_removed_ids.add (removed_id);

				// Also remove it from all_absent_ids because it's an invariant violation to be in both

				all_absent_ids.remove (removed_id);
			}
		}

		return;
	}


	// Create a new empty assignment, and make it the successor to the given predecessor.
	// Parameters:
	//  predecessor = An assignment that is NOT in this list.  It must have a timeline ID.
	//                It must not currently have a successor.
	// Returns the newly-created assignment.

	public AliasAssignment make_empty_successor (AliasAssignment predecessor) {

		// Create the successor

		AliasAssignment successor = new AliasAssignment();

		// Add the successor to the list

		add_assignment (successor);

		// Match it

		match_successor (predecessor, successor);

		// Return it

		return successor;
	}


	// Get the unmatched assignment for a given Comcat ID.
	// Returns null if none found.
	// Note: Returns null if the assignment exists but is already matched.

	public AliasAssignment get_unmatched_assignment_for_comcat_id (String comcat_id) {
		AliasAssignment aa = all_comcat_ids.get (comcat_id);
		if (aa != null) {
			if (aa.has_match()) {
				aa = null;
			}
		}
		return aa;
	}


	// Get the unmatched assignment for a given primary ID.
	// Returns null if none found.
	// Note: Returns null if the assignment exists but is already matched.

	public AliasAssignment get_unmatched_assignment_for_primary_id (String primary_id) {
		AliasAssignment aa = all_comcat_ids.get (primary_id);
		if (aa != null) {
			if (!( aa.is_unmatched() && primary_id.equals (aa.get_primary_id()) )) {		// equals is always false if its argument is null
				aa = null;
			}
		}
		return aa;
	}


	// Return true if there are no Comcat IDs in common between this list and the given assignment.

	public boolean is_disjoint_assignment (AliasAssignment aa) {
		for (Iterator<String> it = aa.get_comcat_id_iterator(); it.hasNext(); ) {
			if (all_comcat_ids.get (it.next()) != null) {
				return false;
			}
		}
		return true;
	}


	// Set the timeline ID for the given assignment.
	// Parameters:
	//  aa = An assignment in this list which does not currently have a timeline ID.
	//  id = The timeline ID to assign.

	public boolean set_timeline_id (AliasAssignment aa, String id) {
		if (!( id != null && aa != null )) {
			throw new IllegalArgumentException ("AliasAssignmentList.set_timeline_id: Null argument supplied");
		}
		if (aa.has_timeline_id()) {
			throw new IllegalArgumentException ("AliasAssignmentList.set_timeline_id: Alias assignment already has a timeline ID");
		}

		aa.set_timeline_id (id);

		if (all_timeline_ids.put (id, aa) != null) {
			throw new IllegalArgumentException ("AliasAssignmentList.set_timeline_id: Duplicate timeline ID: " + id);
		}

		return true;
	}


	//----- Testing -----


	// Check invariant, throw exception if violated.

	public void check_invariant () {

		// Make a set containing all assignments

		HashSet<AliasAssignment> all_assignments = new HashSet<AliasAssignment>();

		// Make a set containing the union of all removed IDs

		HashSet<String> union_removed_ids = new HashSet<String>();

		// Loop over assignments

		for (AliasAssignment aa : assignments) {

			// Add to the set of all assignments

			if (!( all_assignments.add (aa) )) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Duplicate assignment");
			}

			// Check its own invariant

			aa.check_invariant();

			// If it has a timeline ID, check that it is in the map
			// (this also detects duplicate timeline IDs)

			if (aa.has_timeline_id()) {
				if (get_assignment_for_timeline_id (aa.get_timeline_id()) != aa) {
					throw new RuntimeException("AliasAssignmentList.check_invariant: Unmapped timeline ID: " + aa.get_timeline_id());
				}
			}

			// Check that each of its Comcat IDs is in the map
			// (this also detects duplicate Comcat IDs)

			for (Iterator<String> it = aa.get_comcat_id_iterator(); it.hasNext(); ) {
				String comcat_id = it.next();
				if (get_assignment_for_comcat_id (comcat_id) != aa) {
					throw new RuntimeException("AliasAssignmentList.check_invariant: Unmapped Comcat ID: " + comcat_id);
				}
			}

			// Check that each of its removed IDs is in one of the maps

			for (Iterator<String> it = aa.get_removed_id_iterator(); it.hasNext(); ) {
				String removed_id = it.next();
				if (!( is_known_id (removed_id) )) {
					throw new RuntimeException("AliasAssignmentList.check_invariant: Unmapped removed ID: " + removed_id);
				}
				union_removed_ids.add (removed_id);
			}
		}

		// Loop over timeline IDs

		for (String timeline_id : all_timeline_ids.keySet()) {
		
			// Check that it refers to one of our assignments

			AliasAssignment aa = get_assignment_for_timeline_id (timeline_id);
			if (!( all_assignments.contains (aa) )) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Unknown assignment for timeline ID: " + timeline_id);
			}

			// Check the timeline ID matches

			if (!( timeline_id.equals (aa.get_timeline_id()) )) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Mismatch for timeline ID: " + timeline_id);
			}
		}

		// Loop over Comcat IDs

		for (String comcat_id : all_comcat_ids.keySet()) {
		
			// Check that it refers to one of our assignments

			AliasAssignment aa = get_assignment_for_comcat_id (comcat_id);
			if (!( all_assignments.contains (aa) )) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Unknown assignment for Comcat ID: " + comcat_id);
			}

			// Check the assignment contains the Comcat ID

			if (!( aa.contains_comcat_id (comcat_id) )) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Mismatch for Comcat ID: " + comcat_id);
			}
		}

		// Loop over removed IDs

		for (String removed_id : all_removed_ids) {
		
			// Check that it is not also a Comcat ID

			if (get_assignment_for_comcat_id (removed_id) != null) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Removed ID is also a Comcat ID: " + removed_id);
			}

			// Check that it is a removed ID for at least one assignment

			if (!( union_removed_ids.contains (removed_id) )) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: No assignment for removed ID: " + removed_id);
			}
		}

		// Loop over absent IDs

		for (String absent_id : all_absent_ids) {
		
			// Check that it is not also a Comcat ID or removed ID

			if (is_known_id (absent_id)) {
				throw new RuntimeException("AliasAssignmentList.check_invariant: Absent ID is also a known ID: " + absent_id);
			}
		}
	
		return;
	}


	// Check invariant, throw exception if violated.
	// These additional tests should be satisfied by database entries.

	public void check_invariant_db () {

		// Maps IDs to a set of all assignments referencing the ID

		HashMap<String, HashSet<AliasAssignment>> id_to_aa_map = new HashMap<String, HashSet<AliasAssignment>>();
		
		for (String comcat_id : all_comcat_ids.keySet()) {
			id_to_aa_map.put (comcat_id, new HashSet<AliasAssignment>());
		}
		for (String removed_id : all_removed_ids) {
			id_to_aa_map.put (removed_id, new HashSet<AliasAssignment>());
		}

		// Make a set containing all assignments found so far

		HashSet<AliasAssignment> found_assignments = new HashSet<AliasAssignment>();

		// Make a queue of assignments to be processed

		ArrayDeque<AliasAssignment> aa_queue = new ArrayDeque<AliasAssignment>();

		// Loop over assignments

		for (AliasAssignment aa : assignments) {

			// Check its own invariant

			aa.check_invariant_db();

			// Add each Comcat ID to the map

			for (Iterator<String> it = aa.get_comcat_id_iterator(); it.hasNext(); ) {
				String comcat_id = it.next();
				id_to_aa_map.get(comcat_id).add(aa);
			}

			// Add each removed ID to the map

			for (Iterator<String> it = aa.get_removed_id_iterator(); it.hasNext(); ) {
				String removed_id = it.next();
				id_to_aa_map.get(removed_id).add(aa);
			}
		}

		// Put the first assignment on the queue

		if (get_assignment_count() > 0) {
			aa_queue.add (assignments.get(0));
		}

		// Loop until queue is empty

		for (AliasAssignment aa = aa_queue.poll(); aa != null; aa = aa_queue.poll()) {

			// If not previously found ...

			if (found_assignments.add (aa)) {

				// Add to queue all assignments that share one of our Comcat IDs

				for (Iterator<String> it = aa.get_comcat_id_iterator(); it.hasNext(); ) {
					String comcat_id = it.next();
					for (AliasAssignment aa2 : id_to_aa_map.get(comcat_id)) {
						aa_queue.add (aa2);
					}
				}

				// Add to queue all assignments that share one of our removed IDs

				for (Iterator<String> it = aa.get_removed_id_iterator(); it.hasNext(); ) {
					String removed_id = it.next();
					for (AliasAssignment aa2 : id_to_aa_map.get(removed_id)) {
						aa_queue.add (aa2);
					}
				}
			}
		}

		// If we didn't find all assignments, then this is not a family

		if (found_assignments.size() != get_assignment_count()) {
			throw new RuntimeException("AliasAssignmentList.check_invariant_db: Assignment list does not form a family");
		}

		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 35001;

	private static final String M_VERSION_NAME = "AliasAssignmentList";

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Only marshal objects that satisfy the database invariant

		check_invariant();
		check_invariant_db();

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		int n = get_assignment_count();
		writer.marshalArrayBegin ("assignments", n);
		for (AliasAssignment aa : assignments) {
			writer.marshalMapBegin (null);
			writer.marshalString ("timeline_id", aa.get_timeline_id());
			writer.marshalStringArray ("comcat_ids", aa.get_comcat_ids_as_array());
			writer.marshalStringArray ("removed_ids", aa.get_removed_ids_as_array());
			writer.marshalMapEnd ();
		}
		writer.marshalArrayEnd ();
		writer.marshalLong ("family_time", family_time);
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {

		// Clear the object

		clear();
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		int n = reader.unmarshalArrayBegin ("assignments");
		for (int i = 0; i < n; ++i) {
			AliasAssignment aa = new AliasAssignment();
			reader.unmarshalMapBegin (null);
			aa.set_timeline_id (reader.unmarshalString ("timeline_id"));
			aa.set_comcat_ids_from_array (reader.unmarshalStringArray ("comcat_ids"));
			aa.set_removed_ids_from_array (reader.unmarshalStringArray ("removed_ids"));
			reader.unmarshalMapEnd ();
			add_assignment (aa);
		}
		reader.unmarshalArrayEnd ();
		family_time = reader.unmarshalLong ("family_time");

		// Only marshal objects that satisfy the database invariant

		check_invariant();
		check_invariant_db();

		return;
	}

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public AliasAssignmentList unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}



}

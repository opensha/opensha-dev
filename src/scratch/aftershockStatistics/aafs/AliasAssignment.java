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


/**
 * Assignment of Comcat IDs to a timeline.
 * Author: Michael Barall 06/30/2018.
 */
public class AliasAssignment {


	//----- Values -----

	// The timeline ID.
	// This may be null if no timeline ID has been assigned.

	private String timeline_id;

	// The Comcat IDs that are currently assigned to this timeline.

	private HashSet<String> comcat_ids;

	// The primary Comcat ID for this timeline.
	// If comcat_ids is non-empty and a primary ID has been assigned, then this is a string that also appears in comcat_ids.
	// If comcat_ids is empty, or if no primary ID has been assigned, then this is null.

	private String primary_id;

	// The Comcat IDs that were previously assigned to this timeline, but have been removed.
	// Note that comcat_ids and removed_ids are disjoint sets.

	private HashSet<String> removed_ids;

	// This is used when matching assignments between two lists of assignments.
	// If match_assignment is non-null, then it is the matching assignment in the other list.
	// If match_assignment is null, then this assignment is not matched yet.
	// A matching assignment may be either a predecessor or a successor.

	private AliasAssignment match_assignment;


	//----- Construction -----

	// The default constructor creates an empty assignment.

	public AliasAssignment () {
		this.timeline_id = null;
		this.comcat_ids = new HashSet<String>();
		this.primary_id = null;
		this.removed_ids = new HashSet<String>();
		this.match_assignment = null;
	}

	// Clear the assignment.

	public void clear () {
		this.timeline_id = null;
		this.comcat_ids = new HashSet<String>();
		this.primary_id = null;
		this.removed_ids = new HashSet<String>();
		this.match_assignment = null;
	}


	//----- Access -----

	// Get the timeline ID, or null if none has been assigned.

	public String get_timeline_id () {
		return timeline_id;
	}


	// Return true if a timeline ID has been assigned.

	public boolean has_timeline_id () {
		return timeline_id != null;
	}


	// Set the timeline ID.

	public void set_timeline_id (String timeline_id) {
		this.timeline_id = timeline_id;
	}


	// Get the number of Comcat IDs.

	public int get_comcat_id_count () {
		return comcat_ids.size();
	}


	// Get an iterator over Comcat IDs.

	public Iterator<String> get_comcat_id_iterator () {
		return comcat_ids.iterator();
	}


	// Return true if there is at least one Comcat ID.

	public boolean has_comcat_id () {
		return comcat_ids.size() > 0;
	}


	// Return true if the given ID is in the list of Comcat IDs.

	public boolean contains_comcat_id (String id) {
		return comcat_ids.contains(id);
	}


	// Get the primary ID, or null if none has been assigned.

	public String get_primary_id () {
		return primary_id;
	}


	// Return true if a primary ID has been assigned.

	public boolean has_primary_id () {
		return primary_id != null;
	}


	// Get the number of removed IDs.

	public int get_removed_id_count () {
		return removed_ids.size();
	}


	// Get an iterator over removed IDs.

	public Iterator<String> get_removed_id_iterator () {
		return removed_ids.iterator();
	}


	// Return true if there is at least one removed ID.

	public boolean has_removed_id () {
		return removed_ids.size() > 0;
	}


	// Return true if the given ID is in the list of removed IDs.

	public boolean contains_removed_id (String id) {
		return removed_ids.contains(id);
	}


	// Return true if the given ID is in the list of Comcat IDs or removed IDs.

	public boolean is_known_id (String id) {
		return (comcat_ids.contains(id) || removed_ids.contains(id));
	}


//	// Set the Comcat IDs and primary ID from an array.
//	// The first element of the array is the primary ID.
//	// The array must be non-null, have at least one element,
//	// have all element non-null, and have no duplicate elements.
//
//	public void set_comcat_ids (String[] idlist) {
//
//		// Check there is at least one array element
//
//		if (!( idlist != null && idlist.length > 0 )) {
//			throw new IllegalArgumentException("AliasAssignment.set_comcat_ids: No IDs supplied");
//		}
//
//		// Add the array elements to the set of Comcat IDs 
//
//		for (String id: idlist) {
//			if (!( id != null )) {
//				throw new IllegalArgumentException("AliasAssignment.set_comcat_ids: Null ID supplied");
//			}
//			if (!( comcat_ids.add(id) )) {
//				throw new IllegalArgumentException("AliasAssignment.set_comcat_ids: Duplicate ID supplied : " + id);
//			}
//		}
//
//		// Save the primary ID
//
//		primary_id = idlist[0];
//		return;
//	}


	// Add the given ID to the set of Comcat IDs.
	// If f_primary is true and no primary ID is currently defined, then the ID also becomes the primary ID.
	// If the ID is already in comcat_ids or removed_ids, an exception is thrown.
	// Note this routine never changes an existing primary ID.

	public void add_comcat_id (String id, boolean f_primary) {

		// Check arguments

		if (!( id != null )) {
			throw new IllegalArgumentException("AliasAssignment.add_comcat_id: Null ID supplied");
		}

		// Check if ID is in removed_ids

		if (removed_ids.contains(id)) {
			throw new IllegalArgumentException("AliasAssignment.add_comcat_id: Duplicate ID supplied : " + id);
		}

		// Add to comcat_ids, and check if it's a duplicate

		if (!( comcat_ids.add(id) )) {
			throw new IllegalArgumentException("AliasAssignment.add_comcat_id: Duplicate ID supplied : " + id);
		}

		// If we want to be primary ID and there is no primary ID, make it the primary ID

		if (f_primary && primary_id == null) {
			primary_id = id;
		}

		return;
	}


	// Add the given ID to the set of Comcat IDs.
	// If f_primary is true and no primary ID is currently defined, then the ID also becomes the primary ID.
	// This function permits the ID to be a member of comcat_ids or removed_ids.  The ID is
	// removed from removed_ids (if present) and added to comcat_ids (if not present).
	// Note this routine never changes an existing primary ID.

	public void add_comcat_id_dup_ok (String id, boolean f_primary) {

		// Check arguments

		if (!( id != null )) {
			throw new IllegalArgumentException("AliasAssignment.add_comcat_id_dup_ok: Null ID supplied");
		}

		// Remove ID from removed_ids, if it is there

		removed_ids.remove(id);

		// Add to comcat_ids, if it is not there

		comcat_ids.add(id);

		// If we want to be primary ID and there is no primary ID, make it the primary ID

		if (f_primary && primary_id == null) {
			primary_id = id;
		}

		return;
	}


	// Set the Comcat IDs and primary ID from an array range.
	// The range extends from index lo (inclusive) to index hi (exclusive).
	// If lo == hi then the range is empty, in which case nothing is added.
	// If the range is non-empty, then the first element of the range is the primary ID,
	// provided that there is currently no primary ID.
	// The array must be non-null, must satisfy 0 <= lo <= hi <= idlist.length,
	// and all elements in the range must be non-null.
	// An exception is thrown if any ID is already present in comcat_ids or removed_ids.
	// If there are already Comcat IDs defined, these are added to the set.

	public void set_comcat_ids_from_array (String[] idlist, int lo, int hi) {

		// Check arguments

		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignment.set_comcat_ids_from_array: No IDs supplied");
		}

		if (!( 0 <= lo && lo <= hi && hi <= idlist.length )) {
			throw new IllegalArgumentException("AliasAssignment.set_comcat_ids_from_array: Invalid range : lo = " + lo + ", hi = " + hi + ", length = " + idlist.length);
		}

		// Add the array elements to the set of Comcat IDs 

		for (int i = lo; i < hi; ++i) {
			add_comcat_id (idlist[i], true);
		}

		return;
	}


	// Set the Comcat IDs and primary ID from an array.
	// If the array is empty then nothing is added.
	// If the array is non-empty, then the first element of the array is the primary ID,
	// provided that there is currently no primary ID.
	// The array must be non-null, and all elements must be non-null.
	// An exception is thrown if any ID is already present in comcat_ids or removed_ids.
	// If there are already Comcat IDs defined, these are added to the set.

	public void set_comcat_ids_from_array (String[] idlist) {
		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignment.set_comcat_ids_from_array: No IDs supplied");
		}

		set_comcat_ids_from_array (idlist, 0, idlist.length);
		return;
	}


	// Get the Comcat IDs and primary ID into an array range.
	// The range extends from index lo (inclusive) to index hi (exclusive),
	// where hi = lo + get_comcat_id_count().
	// If lo == hi then the range is empty, in which case nothing is retrieved.
	// If the range is non-empty, then the first element of the range is the primary ID.
	// The array must be non-null, and must satisfy 0 <= lo <= hi <= idlist.length.
	// The return value is hi.

	public int get_comcat_ids_as_array (String[] idlist, int lo) {

		// Check arguments

		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignment.get_comcat_ids_as_array: No array supplied");
		}

		int hi = lo + comcat_ids.size();

		if (!( 0 <= lo && lo <= hi && hi <= idlist.length )) {
			throw new IllegalArgumentException("AliasAssignment.get_comcat_ids_as_array: Invalid range : lo = " + lo + ", hi = " + hi + ", length = " + idlist.length);
		}

		// Insert the primary ID, if any

		int i = lo;
		if (primary_id != null) {
			idlist[i++] = primary_id;
		}

		// Insert the Comcat IDs

		for (String id : comcat_ids) {
			if (!( id.equals (primary_id) )) {		// equals is always false if primary_id is null
				idlist[i++] = id;
			}
		}

		return hi;
	}


	// Get the Comcat IDs and primary ID into an array.
	// The length of the returned array equals get_comcat_id_count().
	// If the array is non-empty, then the first element of the array is the primary ID.

	public String[] get_comcat_ids_as_array () {
		String[] idlist = new String[comcat_ids.size()];
		get_comcat_ids_as_array (idlist, 0);
		return idlist;
	}


	// Add the given ID to the set of removed IDs.
	// If the ID is already in comcat_ids or removed_ids, an exception is thrown.

	public void add_removed_id (String id) {

		// Check arguments

		if (!( id != null )) {
			throw new IllegalArgumentException("AliasAssignment.add_removed_id: Null ID supplied");
		}

		// If already in comcat_ids, it's a duplicate

		if (comcat_ids.contains(id)) {
			throw new IllegalArgumentException("AliasAssignment.add_removed_id: Duplicate ID supplied : " + id);
		}

		// Add to removed_ids, and check if it's a duplicate

		if (!( removed_ids.add(id) )) {
			throw new IllegalArgumentException("AliasAssignment.add_removed_id: Duplicate ID supplied : " + id);
		}

		return;
	}


	// Add the given ID to the set of removed IDs.
	// This function permits the ID to be a member of comcat_ids or removed_ids.  In either
	// case, the function performs no operation.

	public void add_removed_id_dup_ok (String id) {

		// Check arguments

		if (!( id != null )) {
			throw new IllegalArgumentException("AliasAssignment.add_removed_id_dup_ok: Null ID supplied");
		}

		// If already in comcat_ids, do nothing

		if (comcat_ids.contains(id)) {
			return;
		}

		// Add to removed_ids, if it's not there

		removed_ids.add(id);

		return;
	}


	// Set the removed IDs from an array range.
	// The range extends from index lo (inclusive) to index hi (exclusive).
	// If lo == hi then the range is empty, in which case nothing is added.
	// The array must be non-null, must satisfy 0 <= lo <= hi <= idlist.length,
	// and all elements in the range must be non-null.
	// An exception is thrown if any ID is already present in comcat_ids or removed_ids.
	// If there are already removed IDs defined, these are added to the set.

	public void set_removed_ids_from_array (String[] idlist, int lo, int hi) {

		// Check arguments

		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignment.set_removed_ids_from_array: No IDs supplied");
		}

		if (!( 0 <= lo && lo <= hi && hi <= idlist.length )) {
			throw new IllegalArgumentException("AliasAssignment.set_removed_ids_from_array: Invalid range : lo = " + lo + ", hi = " + hi + ", length = " + idlist.length);
		}

		// Add the array elements to the set of removed IDs 

		for (int i = lo; i < hi; ++i) {
			add_removed_id (idlist[i]);
		}

		return;
	}


	// Set the removed IDs from an array.
	// If the array is empty then nothing is added.
	// The array must be non-null, and all elements must be non-null.
	// An exception is thrown if any ID is already present in comcat_ids or removed_ids.
	// If there are already removed IDs defined, these are added to the set.

	public void set_removed_ids_from_array (String[] idlist) {
		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignment.set_removed_ids_from_array: No IDs supplied");
		}

		set_removed_ids_from_array (idlist, 0, idlist.length);
		return;
	}


	// Get the removed IDs into an array range.
	// The range extends from index lo (inclusive) to index hi (exclusive),
	// where hi = lo + get_removed_id_count().
	// If lo == hi then the range is empty, in which case nothing is retrieved.
	// The array must be non-null, and must satisfy 0 <= lo <= hi <= idlist.length.
	// The return value is hi.

	public int get_removed_ids_as_array (String[] idlist, int lo) {

		// Check arguments

		if (!( idlist != null )) {
			throw new IllegalArgumentException("AliasAssignment.get_removed_ids_as_array: No array supplied");
		}

		int hi = lo + removed_ids.size();

		if (!( 0 <= lo && lo <= hi && hi <= idlist.length )) {
			throw new IllegalArgumentException("AliasAssignment.get_removed_ids_as_array: Invalid range : lo = " + lo + ", hi = " + hi + ", length = " + idlist.length);
		}

		// Insert the removed IDs

		int i = lo;
		for (String id : removed_ids) {
			idlist[i++] = id;
		}

		return hi;
	}


	// Get the removed IDs into an array.
	// The length of the returned array equals get_removed_id_count().

	public String[] get_removed_ids_as_array () {
		String[] idlist = new String[removed_ids.size()];
		get_removed_ids_as_array (idlist, 0);
		return idlist;
	}


	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "AliasAssignment:\n"
			+ "\ttimeline_id: " + ((timeline_id == null) ? ("null") : (timeline_id)) + "\n"
			+ "\tprimary_id: " + ((primary_id == null) ? ("null") : (primary_id)) + "\n"
			+ "\tcomcat_ids: " + Arrays.toString (get_comcat_ids_as_array()) + "\n"
			+ "\tremoved_ids: " + Arrays.toString (get_removed_ids_as_array());
		return str;
	}


	// toString - Convert to string.
	// This version inserts an index number into the string.

	public String toString (int index) {
		String str = "AliasAssignment[" + index + "]:\n"
			+ "\ttimeline_id: " + ((timeline_id == null) ? ("null") : (timeline_id)) + "\n"
			+ "\tprimary_id: " + ((primary_id == null) ? ("null") : (primary_id)) + "\n"
			+ "\tcomcat_ids: " + Arrays.toString (get_comcat_ids_as_array()) + "\n"
			+ "\tremoved_ids: " + Arrays.toString (get_removed_ids_as_array());
		return str;
	}


	// Make a string representation on a single line.

	public String one_line_string () {

		// Number of Comcat and removed IDs

		int cc_count = get_comcat_id_count();
		int rm_count = get_removed_id_count();

		int total_count = cc_count + rm_count;

		// An array large enough for both

		String[] all_ids = new String[total_count];

		// Fill with IDs, first the Comcat and then the removed IDs

		get_comcat_ids_as_array (all_ids, 0);
		get_removed_ids_as_array (all_ids, cc_count);

		// Prefix each removed ID with a minus sign

		for (int i = cc_count; i < total_count; ++i) {
			all_ids[i] = "-" + all_ids[i];
		}

		// Construct the string

		return ((timeline_id == null) ? ("null") : (timeline_id))
				+ " = " + Arrays.toString (all_ids);
	}


	//----- Matching -----


	// Test if another assignment has the same primary and Comcat IDs as this assignment.

	public boolean is_same_primary_comcat_ids (AliasAssignment other) {
		if (!( (primary_id == null) ? (other.primary_id == null) : (primary_id.equals(other.primary_id)) )) {
			return false;
		}
		if (!( comcat_ids.size() == other.comcat_ids.size() )) {
			return false;
		}
		if (!( comcat_ids.containsAll (other.comcat_ids) )) {
			return false;
		}
		return true;
	}


	// Add removed IDs from a predecessor assignment.
	// All IDs in other.comcat_ids and other.removed_ids are added to our removed_ids.
	// If any ID is already in our comcat_ids or removed_ids, then that ID is ignored.

	public void add_removed_ids_from_predecessor (AliasAssignment predecessor) {
		if (!( predecessor != null )) {
			throw new IllegalArgumentException("AliasAssignment.add_removed_ids_from_predecessor: No predecessor assignment supplied");
		}

		for (String id : predecessor.comcat_ids) {
			add_removed_id_dup_ok (id);
		}

		for (String id : predecessor.removed_ids) {
			add_removed_id_dup_ok (id);
		}

		return;
	}


	// Get the matching assignment, or null if none if not matched yet.

	public AliasAssignment get_match_assignment () {
		return match_assignment;
	}


	// Set the matching assignment.

	public void set_match_assignment (AliasAssignment assignment) {

		// Check arguments

		if (!( assignment != null )) {
			throw new IllegalArgumentException ("AliasAssignment.get_match_assignment: Null assignment supplied");
		}

		// Error if already matched

		if (match_assignment != null) {
			throw new IllegalStateException ("AliasAssignment.get_match_assignment: Assignment is already matched");
		}

		// Set the match

		match_assignment = assignment;
		return;
	}


	// Return true if a matching assignment has been set.

	public boolean has_match () {
		return match_assignment != null;
	}


	// Return true if a matching assignment has NOT been set.

	public boolean is_unmatched () {
		return match_assignment == null;
	}


	// Compare the timeline IDs of this assignment and another assignment.
	// Returns negative, zero, or positive depending on whether this assignments's timeline ID
	// is less than, equal to, or greater than the other assignments's timeline ID.
	// Throws an exception if either assignment does not have a timeline ID.

	public int compare_timeline_id_to (AliasAssignment other) {
		if (!( this.timeline_id != null && other.timeline_id != null )) {
			throw new IllegalStateException ("AliasAssignment.compare_timeline_id_to: Encountered assignment with no timeline ID");
		}
		return this.timeline_id.compareTo (other.timeline_id);
	}


	// Compare the primary IDs of this assignment and another assignment.
	// Returns negative, zero, or positive depending on whether this assignments's primary ID
	// is less than, equal to, or greater than the other assignments's primary ID.
	// Throws an exception if either assignment does not have a primary ID.

	public int compare_primary_id_to (AliasAssignment other) {
		if (!( this.primary_id != null && other.primary_id != null )) {
			throw new IllegalStateException ("AliasAssignment.compare_primary_id_to: Encountered assignment with no primary ID");
		}
		return this.primary_id.compareTo (other.primary_id);
	}


	//----- Testing -----


	// Check invariant, throw exception if violated.

	public void check_invariant () {

		// Check that the primary ID, if defined, is a Comcat ID

		if (has_primary_id()) {
			if (!( contains_comcat_id (get_primary_id()) )) {
				throw new RuntimeException("AliasAssignment.check_invariant: Primary ID is not a Comcat ID: " + get_primary_id());
			}
		}

		// Check that removed_ids and comcat_ids are disjoint sets

		for (String id : removed_ids) {
			if (contains_comcat_id (id)) {
				throw new RuntimeException("AliasAssignment.check_invariant: Removed ID is also a Comcat ID: " + id);
			}
		}
	
		return;
	}


	// Check invariant, throw exception if violated.
	// These additional tests should be satisfied by database entries.

	public void check_invariant_db () {

		// Check that there is a timeline ID

		if (!( has_timeline_id() )) {
			throw new RuntimeException("AliasAssignment.check_invariant_db: Timeline ID is not defined");
		}

		// Check there is at least one known ID

		if (!( has_comcat_id() || has_removed_id() )) {
			throw new RuntimeException("AliasAssignment.check_invariant_db: No known IDs");
		}

		// Check that there is a primary ID, if there are any Comcat IDs

		if (has_comcat_id()) {
			if (!( has_primary_id() )) {
				throw new RuntimeException("AliasAssignment.check_invariant_db: Primary ID is not defined, but there are Comcat IDs");
			}
		}
	
		return;
	}


}

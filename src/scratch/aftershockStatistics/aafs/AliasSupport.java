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

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;
import scratch.aftershockStatistics.aafs.entity.AliasFamily;

import scratch.aftershockStatistics.util.EventNotFoundException;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.SimpleUtils;

import scratch.aftershockStatistics.CompactEqkRupList;
import scratch.aftershockStatistics.ComcatException;
import scratch.aftershockStatistics.ComcatConflictException;
import scratch.aftershockStatistics.ComcatRemovedException;

/**
 * Support functions for aliases.
 * Author: Michael Barall 07/07/2018.
 *
 * This class contains functions for maintaining the mapping from Comcat IDs to timeline IDs.
 */
public class AliasSupport extends ServerComponent {




	//----- Alias configuration -----

	// True to enable verbose mode.

	private boolean alias_verbose = true;




	//----- Timeline ID -----
	//
	// Each timeline is assigned an ID that is distinct from any possible Comcat ID.
	// That ID is used as the "event ID" for the timeline.
	// This allows a timeline to retain a constant ID even if the Comcat IDs change.
	//
	// Ordinarily, the timeline ID is the primary Comcat ID with a dot prepended.
	// So Comcat event "us2000abcd" would have timeline ID ".us2000abcd".  If for some
	// reason ".us2000abcd" is not available (which might happen if Comcat IDs have
	// been changing) then we use ".us2000abcd.2" or ".us2000abcd.3" etc.  For a given
	// timeline, the timeline ID remains constant even if Comcat changes the primary ID.
	//
	// The PDL code is obtained by removing the initial dot from the timeline ID.
	// So in most cases, the PDL code is the Comcat primary ID.




	// The prefix used for timeline IDs.

	public static final String TLID_PREFIX = ".";

	// The separator used when appending numbers to the end of a timeline ID.

	public static final String TLID_SEPARATOR = ".";

	// Maximum suffix allowed in a timeline ID.

	public static final int TLID_MAX_SUFFIX = 999999;




	// Return true if the given ID is a timeline id.

	public boolean is_timeline_id (String id) {
		return id.startsWith (TLID_PREFIX);
	}



	// Given a timeline ID, return the PDL code.

	public String timeline_id_to_pdl_code (String timeline_id) {
		if (!( timeline_id.startsWith (TLID_PREFIX) )) {
			throw new IllegalArgumentException("AliasSupport.timeline_id_to_pdl_code: Invalid timeline ID supplied: " + timeline_id);
		}

		return timeline_id.substring (TLID_PREFIX.length());
	}




	// Create a new timeline ID, given the Comcat primary ID.
	// This function checks that the new timeline ID does not appear in the database.
	// Note: This function must be designed so that calls with two different Comcat IDs
	// can never return the same timeline ID regardless of database contents.

	public String create_timeline_id (String comcat_id) {
	
		// Loop over possible suffixes

		for (int suffix = 1; suffix <= TLID_MAX_SUFFIX; ++suffix) {
		
			// Construct the possible timeline ID

			String timeline_id;
			if (suffix == 1) {
				timeline_id = TLID_PREFIX + comcat_id;
			} else {
				timeline_id = TLID_PREFIX + comcat_id + TLID_SEPARATOR + suffix;
			}

			// Check if the timeline ID appears in the database

			AliasFamily alfam = AliasFamily.get_recent_alias_family (0L, 0L, timeline_id, null, null);

			// If not found, then this timeline ID is OK

			if (alfam == null) {

				if (alias_verbose) {
					System.out.println ("Alias: Created new timeline ID = " + timeline_id + ", for Comcat ID = " + comcat_id);
				}

				return timeline_id;
			}
		}

		// Exhausted possible suffixes

		throw new DBCorruptException ("AliasSupport.create_timeline_id - Exhaused all possible timeline IDs");
	}




	//----- Alias tracking -----




	// Obtain assignment lists from the database, and merge them.
	// Parameters:
	//  assignments = Existing assignment list, to receive merge.
	//  in_queue = Input queue of Comcat IDs to process.
	//  out_queue = Output queue of Comcat IDs that require further processing.
	// The return value is the number of merges performed.
	// This function retrieves Comcat IDs from in_queue until the queue is empty.
	// For each ID that is not already known to assignments, fetch the corresponding
	// assigment list from the database, and merge it into assignments.  Also, all
	// Comcat IDs (both current and removed) in the fetched list are added to out_queue.
	// An exception thrown from this function should be considered a database error.

	private int merge_from_database (AliasAssignmentList assignments, Queue<String> in_queue, Queue<String> out_queue) {
	
		int merge_count = 0;
		int id_count = 0;

		// Loop until input queue is empty

		for (String id = in_queue.poll(); id != null; id = in_queue.poll()) {

			++id_count;

			// If the ID is not known ...

			if (!( assignments.is_known_or_absent_id (id) )) {

				// Fetch database entry for this ID

				String[] ids = new String[1];
				ids[0] = id;
				AliasFamily alfam = AliasFamily.get_recent_alias_family (0L, 0L, null, ids, null);

				// If we got one ...

				if (alfam != null) {

					// Get the new assignment list

					AliasAssignmentList new_assigments = alfam.get_assignments();

					// Add all the new IDs to the output queue
				
					for (Iterator<String> it = new_assigments.get_all_comcat_id_iterator(); it.hasNext(); ) {
						String id2 = it.next();
						if (assignments.contains_absent_id (id2)) {
							throw new DBCorruptException ("AliasSupport.merge_from_database: Encountered absent ID: " + id2);
						}
						out_queue.add (id2);
					}
				
					for (Iterator<String> it = new_assigments.get_all_removed_id_iterator(); it.hasNext(); ) {
						String id2 = it.next();
						if (assignments.contains_absent_id (id2)) {
							throw new DBCorruptException ("AliasSupport.merge_from_database: Encountered absent ID: " + id2);
						}
						out_queue.add (id2);
					}

					// Merge the new assignment list

					assignments.merge_from (new_assigments);
					++merge_count;
				}

				// Otherwise, record it as an absent ID

				else {
					assignments.add_absent_id (id);
				}
			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Merges from database = " + merge_count + ", IDs processed = " + id_count);
		}

		return merge_count;
	}




	// Obtain assignment lists from Comcat, and merge them.
	// Parameters:
	//  assignments = Existing assignment list, to receive merge.
	//  in_queue = Input queue of Comcat IDs to process.
	//  out_queue = Output queue of Comcat IDs that require further processing.
	// The return value is the number of merges performed.
	// This function retrieves Comcat IDs from in_queue until the queue is empty.
	// For each ID that is not already known to assignments, fetch the corresponding
	// assigment list from Comcat, and merge it into assignments.  Also, all
	// Comcat IDs (both current and removed) in the fetched list are added to out_queue.
	// An exception thrown from this function should be considered a Comcat error
	// (which means the operation should be retried later).

	private int merge_from_comcat (AliasAssignmentList assignments, Queue<String> in_queue, Queue<String> out_queue) {
	
		int merge_count = 0;
		int id_count = 0;
	
		// Loop until input queue is empty

		for (String id = in_queue.poll(); id != null; id = in_queue.poll()) {

			++id_count;

			// If the ID is not known ...

			if (!( assignments.is_known_or_absent_id (id) )) {

				// Fetch information from Comcat

				ForecastMainshock fcmain = new ForecastMainshock();
				fcmain.setup_mainshock_poll (id);

				// If we got it ...

				if (fcmain.mainshock_avail) {

					// Get the new assignment, which has no timeline ID or removed IDs

					AliasAssignment new_assigment = new AliasAssignment();
					new_assigment.set_comcat_ids_from_array (fcmain.mainshock_id_list);

					// Add all the new IDs to the output queue
				
					for (Iterator<String> it = new_assigment.get_comcat_id_iterator(); it.hasNext(); ) {
						String id2 = it.next();
						if (assignments.contains_absent_id (id2)) {
							throw new ComcatConflictException ("AliasSupport.merge_from_comcat: Encountered absent ID: " + id2);
						}
						out_queue.add (id2);
					}

					// Add the new assignment

					assignments.add_assignment (new_assigment);
					++merge_count;
				}

				// Otherwise, record it as an absent ID

				else {
					assignments.add_absent_id (id);
				}
			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Merges from Comcat = " + merge_count + ", IDs processed = " + id_count);
		}

		return merge_count;
	}




	// Merge assigment lists from database and Comcat.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  These assignments
	//    must not have timeline IDs or removed IDs.
	// If db_assignments and cc_assignments are both non-empty on entry, then each assignment
	// in cc_assignments must contain at least one Comcat ID that also appears in db_assignments
	// (as either a Comcat ID or a removed ID).
	// This function does two things:
	// - Merge assignments from Comcat into cc_assignments, until cc_assignments contains every
	//   Comcat ID that appears in db_assignments (except those that are unknown to Comcat).
	// - Merge assignments from the database into db_assignments, until db_assignments contains
	//   every Comcat ID that appears in cc_assigments (except those that appear nowhere in the
	//   database).
	// The two merges are performed repeatedly until no further merges are possible.  It is
	// expected that typically this will require about 2 cycles.  The merges preserve the
	// conditions on db_assignments and cc_assignments as described above.
	// The function throws ComcatException if the operation fails due to a problem accessing
	// Comcat (in which case the function should be retried later).  Any other exception
	// likely indicates a problem with the database.

	private void dual_merge_database_comcat (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments) {
	
		// Input queue for Comcat merges

		ArrayDeque<String> cc_in_queue = new ArrayDeque<String>();

		// Initialize it with all the IDs from the database
				
		for (Iterator<String> it = db_assignments.get_all_comcat_id_iterator(); it.hasNext(); ) {
			cc_in_queue.add (it.next());
		}
				
		for (Iterator<String> it = db_assignments.get_all_removed_id_iterator(); it.hasNext(); ) {
			cc_in_queue.add (it.next());
		}
	
		// Input queue for database merges

		ArrayDeque<String> db_in_queue = new ArrayDeque<String>();

		// Initialize it with all the IDs from Comcat
				
		for (Iterator<String> it = cc_assignments.get_all_comcat_id_iterator(); it.hasNext(); ) {
			db_in_queue.add (it.next());
		}
				
		for (Iterator<String> it = cc_assignments.get_all_removed_id_iterator(); it.hasNext(); ) {
			db_in_queue.add (it.next());
		}

		// Deadman counter

		int deadman = 100;		// maximum cycles

		// Loop until both queues are empty

		while (!( cc_in_queue.isEmpty() && db_in_queue.isEmpty() )) {

			// Check the deadman counter

			--deadman;
			if (deadman == 0) {
				throw new DBCorruptException ("AliasSupport.dual_merge_database_comcat: Merge operations not converging");
			}

			// Merge from Comcat

			try {
				merge_from_comcat (cc_assignments, cc_in_queue, db_in_queue);
			}
			catch (ComcatConflictException e) {
				throw new ComcatConflictException ("AliasSupport.dual_merge_database_comcat: Comcat error while merging", e);
			}
			catch (ComcatException e) {
				throw new ComcatException ("AliasSupport.dual_merge_database_comcat: Comcat error while merging", e);
			}
			catch (Exception e) {
				throw new ComcatConflictException ("AliasSupport.dual_merge_database_comcat: Comcat error while merging", e);
			}

			// Merge from database

			try {
				merge_from_database (db_assignments, db_in_queue, cc_in_queue);
			}
			catch (Exception e) {
				throw new DBCorruptException ("AliasSupport.dual_merge_database_comcat: Database error while merging", e);
			}
		}

		return;
	}




	// Match timelines from database and Comcat, method 1.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function attempts to create matches between a predecessor in db_assignments
	// and a successor in cc_assignments.
	// The function looks for pairs (pred, succ) satisfying:
	// - The primary ID of pred is a Comcat ID of succ.
	// - The primary ID of succ is a Comcat ID of pred.
	// These conditions will handle the following simple cases:
	// - Pass-thru: pred and succ have the same primary and Comcat IDs.
	// - Secondary ID change: pred and succ have the same primary ID, but may have different
	//   secondary IDs.
	// - Authoritative change: pred and succ have the same Comcat IDs, but different primary IDs.
	// - Simple merge: succ's Comcat IDs are the union of the Comcat IDs of several possible
	//   predecessors; pred is taken to be the one that supplied succ's primary ID.
	// - Simple split: pred's Comcat IDs are the union of the Comcat IDs of several possible
	//   successors; succ is taken to be the one that received pred's primary ID.
	// Some similar cases involving the addition or deletion of secondary IDs are also handled.

	private int match_timelines_1 (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in cc_assignments
		// (Note that in cc_assignments each assignment is guaranteed to have a primary ID)

		for (Iterator<AliasAssignment> itsucc = cc_assignments.get_assigment_iterator(); itsucc.hasNext(); ) {
			AliasAssignment succ = itsucc.next();
			if (succ.is_unmatched()) {

				// An unmatched element of db_assignments that contains succ's primary ID is a possible predecessor

				AliasAssignment pred = db_assignments.get_unmatched_assignment_for_comcat_id (succ.get_primary_id());
				if (pred != null) {
					
					// [At this point we know that succ's primary ID is one of pred's Comcat IDs]
					// If the predecessor's primary ID (which must exist because it has a Comcat ID) is one of succ's Comcat IDs ...

					if (succ.contains_comcat_id (pred.get_primary_id())) {

						// Match

						cc_assignments.match_successor (pred, succ);
						++match_count;
					}
				}

			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method 1 = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method 2.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	//  f_unique = True to require unique predecessor.  See below.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function attempts to create matches between a predecessor in db_assignments
	// and a successor in cc_assignments.
	// The function looks for pairs (pred, succ) satisfying:
	// - The primary ID of pred is a Comcat ID of succ.
	// - If there is more than one possible predecessor that satisfies the above condition,
	//   then pred is the one with the earliest timeline ID (this ensures a well-defined
	//   result regardless of iteration orders).
	// - If f_unique is true, there must be only one possible predecessor.
	// These conditions will handle the following simple cases:
	// - Superseding authoritative ID: succ contains the Comcat IDs of pred, plus a new ID
	//   which becomes the primary ID.
	// - Superseding split: pred's Comcat IDs are a subset of the union of the Comcat IDs of several
	//   possible successors; succ is taken to be the one that received pred's primary ID.
	// Some similar cases involving the addition or deletion of secondary IDs are also handled.
	// Note: succ is the assignment that would be obtained by simply continuing to query
	// with pred's primary ID.

	private int match_timelines_2 (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments, boolean f_unique) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in cc_assignments

		for (Iterator<AliasAssignment> itsucc = cc_assignments.get_assigment_iterator(); itsucc.hasNext(); ) {
			AliasAssignment succ = itsucc.next();
			if (succ.is_unmatched()) {

				// Set up to search for predecessor

				boolean f_match = false;		// true if match found
				AliasAssignment pred = null;	// the predecessor to use
				HashSet<AliasAssignment> seen = new HashSet<AliasAssignment>();		// assignments already seen

				// Loop over Comcat IDs to find possible predecessors

				for (Iterator<String> it = succ.get_comcat_id_iterator(); it.hasNext(); ) {
					AliasAssignment aa = db_assignments.get_unmatched_assignment_for_primary_id (it.next());
					if (aa != null && seen.add(aa)) {
					
						// [At this point we know that aa's primary ID is one of succ's Comcat IDs]
						// If this is the first possible predecessor, just save it

						if (pred == null) {
							pred = aa;
							f_match = true;
						}

						// Otherwise, this is a new possible predecessor ...

						else {

							// If uniqueness required, no match

							if (f_unique) {
								f_match = false;
							}

							// If this is a new minimum timeline ID, save it

							if (aa.compare_timeline_id_to (pred) < 0) {
								pred = aa;
							}
						}
					}
				}

				// If match is found, match it

				if (f_match) {
					cc_assignments.match_successor (pred, succ);
					++match_count;
				}

			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method 2" + (f_unique ? "U" : "N") + " = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method 3.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	//  f_unique = True to require unique successor.  See below.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function attempts to create matches between a predecessor in db_assignments
	// and a successor in cc_assignments.
	// The function looks for pairs (pred, succ) satisfying:
	// - The primary ID of succ is a Comcat ID of pred.
	// - If there is more than one possible successor that satisfies the above condition,
	//   then succ is the one with the earliest primary ID (this ensures a well-defined
	//   result regardless of iteration orders).
	// - If f_unique is true, there must be only one possible successor.
	// These conditions will handle the following simple cases:
	// - Deleted authoritative ID: The primary ID of pred is deleted, and succ contains the
	//   remaining Comcat IDs of pred.
	// - Superseding merge: succ's Comcat IDs are a subset of the union of the Comcat IDs of several
	//   possible predecessors; pred is taken to be the one that supplied succ's primary ID.
	// Some similar cases involving the addition or deletion of secondary IDs are also handled.

	private int match_timelines_3 (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments, boolean f_unique) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in db_assignments

		for (Iterator<AliasAssignment> itpred = db_assignments.get_assigment_iterator(); itpred.hasNext(); ) {
			AliasAssignment pred = itpred.next();
			if (pred.is_unmatched()) {

				// Set up to search for successor

				boolean f_match = false;		// true if match found
				AliasAssignment succ = null;	// the successor to use
				HashSet<AliasAssignment> seen = new HashSet<AliasAssignment>();		// assignments already seen

				// Loop over Comcat IDs to find possible successors, noting that pred might have no Comcat IDs

				for (Iterator<String> it = pred.get_comcat_id_iterator(); it.hasNext(); ) {
					AliasAssignment aa = cc_assignments.get_unmatched_assignment_for_primary_id (it.next());
					if (aa != null && seen.add(aa)) {
					
						// [At this point we know that aa's primary ID is one of pred's Comcat IDs]
						// If this is the first possible successor, just save it

						if (succ == null) {
							succ = aa;
							f_match = true;
						}

						// Otherwise, this is a new possible successor ...

						else {

							// If uniqueness required, no match

							if (f_unique) {
								f_match = false;
							}

							// If this is a new minimum primary ID, save it

							if (aa.compare_primary_id_to (succ) < 0) {
								succ = aa;
							}
						}
					}
				}

				// If match is found, match it

				if (f_match) {
					cc_assignments.match_successor (pred, succ);
					++match_count;
				}

			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method 3" + (f_unique ? "U" : "N") + " = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method 4.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	//  f_unique = True to require unique predecessor.  See below.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function attempts to create matches between a predecessor in db_assignments
	// and a successor in cc_assignments.
	// The function looks for pairs (pred, succ) satisfying:
	// - There is a Comcat ID in common between pred and succ.
	// - If there is more than one possible predecessor that satisfies the above condition,
	//   then pred is the one with the earliest timeline ID (this ensures a well-defined
	//   result regardless of iteration orders).
	// - If f_unique is true, there must be only one possible predecessor.
	// Note: If there are still unmatched assignments that satisfy the first condition after
	// application of methods 1 thru 3, it is not clear how to pick a "best" pairing as that
	// would indicate the Comcat IDs have been reshuffled.  Methods 4 and 5 repeat methods
	// 2 and 3 but without reference to primary IDs.  It is possible to think of further
	// refinements, but it is not clear that further refinements would be of value.

	private int match_timelines_4 (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments, boolean f_unique) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in cc_assignments

		for (Iterator<AliasAssignment> itsucc = cc_assignments.get_assigment_iterator(); itsucc.hasNext(); ) {
			AliasAssignment succ = itsucc.next();
			if (succ.is_unmatched()) {

				// Set up to search for predecessor

				boolean f_match = false;		// true if match found
				AliasAssignment pred = null;	// the predecessor to use
				HashSet<AliasAssignment> seen = new HashSet<AliasAssignment>();		// assignments already seen

				// Loop over Comcat IDs to find possible predecessors

				for (Iterator<String> it = succ.get_comcat_id_iterator(); it.hasNext(); ) {
					AliasAssignment aa = db_assignments.get_unmatched_assignment_for_comcat_id (it.next());
					if (aa != null && seen.add(aa)) {
					
						// [At this point we know that aa and succ have a Comcat ID in common]
						// If this is the first possible predecessor, just save it

						if (pred == null) {
							pred = aa;
							f_match = true;
						}

						// Otherwise, this is a new possible predecessor ...

						else {

							// If uniqueness required, no match

							if (f_unique) {
								f_match = false;
							}

							// If this is a new minimum timeline ID, save it

							if (aa.compare_timeline_id_to (pred) < 0) {
								pred = aa;
							}
						}
					}
				}

				// If match is found, match it

				if (f_match) {
					cc_assignments.match_successor (pred, succ);
					++match_count;
				}

			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method 4" + (f_unique ? "U" : "N") + " = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method 5.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	//  f_unique = True to require unique successor.  See below.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function attempts to create matches between a predecessor in db_assignments
	// and a successor in cc_assignments.
	// The function looks for pairs (pred, succ) satisfying:
	// - There is a Comcat ID in common between pred and succ.
	// - If there is more than one possible successor that satisfies the above condition,
	//   then succ is the one with the earliest primary ID (this ensures a well-defined
	//   result regardless of iteration orders).
	// - If f_unique is true, there must be only one possible successor.
	// Note: If there are still unmatched assignments that satisfy the first condition after
	// application of methods 1 thru 3, it is not clear how to pick a "best" pairing as that
	// would indicate the Comcat IDs have been reshuffled.  Methods 4 and 5 repeat methods
	// 2 and 3 but without reference to primary IDs.  It is possible to think of further
	// refinements, but it is not clear that further refinements would be of value.

	private int match_timelines_5 (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments, boolean f_unique) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in db_assignments

		for (Iterator<AliasAssignment> itpred = db_assignments.get_assigment_iterator(); itpred.hasNext(); ) {
			AliasAssignment pred = itpred.next();
			if (pred.is_unmatched()) {

				// Set up to search for successor

				boolean f_match = false;		// true if match found
				AliasAssignment succ = null;	// the successor to use
				HashSet<AliasAssignment> seen = new HashSet<AliasAssignment>();		// assignments already seen

				// Loop over Comcat IDs to find possible successors, noting that pred might have no Comcat IDs

				for (Iterator<String> it = pred.get_comcat_id_iterator(); it.hasNext(); ) {
					AliasAssignment aa = cc_assignments.get_unmatched_assignment_for_comcat_id (it.next());
					if (aa != null && seen.add(aa)) {
					
						// [At this point we know that aa and pred have a Comcat ID in common]
						// If this is the first possible successor, just save it

						if (succ == null) {
							succ = aa;
							f_match = true;
						}

						// Otherwise, this is a new possible successor ...

						else {

							// If uniqueness required, no match

							if (f_unique) {
								f_match = false;
							}

							// If this is a new minimum primary ID, save it

							if (aa.compare_primary_id_to (succ) < 0) {
								succ = aa;
							}
						}
					}
				}

				// If match is found, match it

				if (f_match) {
					cc_assignments.match_successor (pred, succ);
					++match_count;
				}

			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method 5" + (f_unique ? "U" : "N") + " = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method revive.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	//  f_unique = True to require unique successor.  See below.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function attempts to create matches between a predecessor in db_assignments
	// and a successor in cc_assignments.
	// The function looks for pairs (pred, succ) satisfying:
	// - There are no Comcat IDs in pred (that is, it is empty).
	// - There is a Comcat ID in succ that is also a removed ID in pred.
	// - There is no Comcat ID in succ that is also a Comcat ID for any assignment in db_assignments.
	// - If there is more than one possible successor that satisfies the above conditions,
	//   then succ is the one with the earliest primary ID (this ensures a well-defined
	//   result regardless of iteration orders).
	// - If f_unique is true, there must be only one possible successor.
	// Note: The purpose of this function is to revive timelines that were previously stopped
	// due to their Comcat IDs being deleted, but which have now re-appeared.

	private int match_timelines_revive (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments, boolean f_unique) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in db_assignments that have no Comcat IDs

		for (Iterator<AliasAssignment> itpred = db_assignments.get_assigment_iterator(); itpred.hasNext(); ) {
			AliasAssignment pred = itpred.next();
			if (pred.is_unmatched() && !(pred.has_comcat_id())) {

				// Set up to search for successor

				boolean f_match = false;		// true if match found
				AliasAssignment succ = null;	// the successor to use
				HashSet<AliasAssignment> seen = new HashSet<AliasAssignment>();		// assignments already seen

				// Loop over removed IDs to find possible successors

				for (Iterator<String> it = pred.get_removed_id_iterator(); it.hasNext(); ) {
					AliasAssignment aa = cc_assignments.get_unmatched_assignment_for_comcat_id (it.next());
					if (aa != null && seen.add(aa) && db_assignments.is_disjoint_assignment(aa)) {
					
						// [At this point we know that aa has a Comcat ID which is a removed ID in pred,
						// and aa has no Comcat ID in common with any assignment in db_assignments]
						// If this is the first possible successor, just save it

						if (succ == null) {
							succ = aa;
							f_match = true;
						}

						// Otherwise, this is a new possible successor ...

						else {

							// If uniqueness required, no match

							if (f_unique) {
								f_match = false;
							}

							// If this is a new minimum primary ID, save it

							if (aa.compare_primary_id_to (succ) < 0) {
								succ = aa;
							}
						}
					}
				}

				// If match is found, match it

				if (f_match) {
					cc_assignments.match_successor (pred, succ);
					++match_count;
				}

			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method R" + (f_unique ? "U" : "N") + " = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method stop.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	// Each assignment in cc_assignments must contain at least one Comcat ID that also
	// appears in db_assignments (as either a Comcat ID or a removed ID).
	// The return value is the number of matches that were made.
	// This function looks for all unmatched assignments in db_assignments.
	// For each of them, create a new empty assignment in cc_assignments to be its successor.

	private int match_timelines_stop (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments) {

		// Number of matches

		int match_count = 0;

		// Loop over unmatched assignments in db_assignments

		for (Iterator<AliasAssignment> itpred = db_assignments.get_assigment_iterator(); itpred.hasNext(); ) {
			AliasAssignment pred = itpred.next();
			if (pred.is_unmatched()) {

				// Create a new empty successor

				cc_assignments.make_empty_successor (pred);
				++match_count;
			}
		}

		if (alias_verbose) {
			System.out.println ("Alias: Matches with method S = " + match_count);
		}

		return match_count;
	}




	// Match timelines from database and Comcat, method main.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  Unmatched assignments
	//    must not have timeline IDs or removed IDs.
	//  new_timelines = List that will receive all new timelines (an assignment that has
	//    Comcat IDs and no predecessor).  These assignments do not have timeline IDs.
	//  revived_timelines = List that will receive all revived timelines (an assignment that has
	//    Comcat IDs, but whose predecessor has no Comcat IDs).
	//  stopped_timelines = List that will receive all stopped timelines (an assignment that has
	//    no Comcat IDs, but whose predecessor has Comcat IDs).
	//  live_timelines = List that will receive all live timelines (an assignment that has
	//    Comcat IDs, and the predecessor also has Comcat IDs).
	//  dead_timelines = List that will receive all dead timelines (an assignment that has
	//    no Comcat IDs, and the predecessor also has no Comcat IDs).
	// Note that all assignments in cc_assignments, except those in new_timelines, have timeline IDs
	// that are already present in the database.  This is required so it is possible to call
	// create_timeline_id() without first writing cc_assignments to the database.
	// Note that every assignment in cc_assignments appears in one of the five lists.

	private void match_timelines_main (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments,
			List<AliasAssignment> new_timelines, List<AliasAssignment> revived_timelines, List<AliasAssignment> stopped_timelines,
			List<AliasAssignment> live_timelines, List<AliasAssignment> dead_timelines) {

		// Try all the methods in turn

		match_timelines_1 (db_assignments, cc_assignments);

		match_timelines_2 (db_assignments, cc_assignments, true);

		match_timelines_3 (db_assignments, cc_assignments, true);

		match_timelines_2 (db_assignments, cc_assignments, false);

		match_timelines_3 (db_assignments, cc_assignments, false);

		match_timelines_4 (db_assignments, cc_assignments, true);

		match_timelines_5 (db_assignments, cc_assignments, true);

		match_timelines_4 (db_assignments, cc_assignments, false);

		match_timelines_5 (db_assignments, cc_assignments, false);

		match_timelines_revive (db_assignments, cc_assignments, true);

		match_timelines_revive (db_assignments, cc_assignments, false);

		match_timelines_stop (db_assignments, cc_assignments);

		// Loop over assignments in cc_assignments

		for (Iterator<AliasAssignment> itsucc = cc_assignments.get_assigment_iterator(); itsucc.hasNext(); ) {
			AliasAssignment aa = itsucc.next();

			// If it is unmatched, then it's new

			if (aa.is_unmatched()) {
				new_timelines.add(aa);
			}

			// Otherwise ...

			else {

				// If it has Comcat IDs ...

				if (aa.has_comcat_id()) {

					// If the predecessor does not have Comcat IDs, then it's revived

					if (!( aa.get_match_assignment().has_comcat_id() )) {
						revived_timelines.add(aa);
					}

					// Otherwise it's live

					else {
						live_timelines.add(aa);
					}
				}

				// Otherwise, it has no Comcat IDs ...

				else {

					// If the predecessor has Comcat IDs, then it's stopped

					if ( aa.get_match_assignment().has_comcat_id() ) {
						stopped_timelines.add(aa);
					}

					// Otherwise it's dead

					else {
						dead_timelines.add(aa);
					}
				}
			}
		}

		return;
	}




	// Delete all notification tasks with the given timeline ID.
	// Note: The currently active task is not deleted, even if it would match.

	public void delete_alias_notification_tasks (String timeline_id) {
		sg.task_sup.delete_all_waiting_tasks (timeline_id, OPCODE_ALIAS_SPLIT, OPCODE_ALIAS_STOP, OPCODE_ALIAS_REVIVE);
		return;
	}




	// Update assigment lists from database and Comcat.
	// Parameters:
	//  db_assignments = List of assignments obtained from the database.  Each assigment
	//    must have a timeline ID, and must have at least one Comcat or removed ID.
	//    Empty assignments (that contain no Comcat IDs) are permitted.  A non-empty
	//    assignment must have a primary ID.
	//  cc_assignments = List of assignments obtained from Comcat.  Each assignment must
	//    contain a primary ID (and therefore at least one Comcat ID).  These assignments
	//    must not have timeline IDs or removed IDs.
	// If db_assignments and cc_assignments are both non-empty on entry, then each assignment
	// in cc_assignments must contain at least one Comcat ID that also appears in db_assignments
	// (as either a Comcat ID or a removed ID).
	// This function does the following:
	// - Merge in assignments from the database and Comcat, as described in dual_merge_database_comcat().
	// - Assign new timeline IDs to the new timelines.
	// - Calculate the family time.
	// - Create tasks for new, revived, and stopped timelines, so timeline status can be updated.
	// - Write cc_assignments to the database as the new family.
	// The function throws ComcatException if the operation fails due to a problem accessing
	// Comcat (in which case the function should be retried later).  Any other exception
	// likely indicates a problem with the database.

	private void update_timelines (AliasAssignmentList db_assignments, AliasAssignmentList cc_assignments) {

		// Merge assignments from database and Comcat

		dual_merge_database_comcat (db_assignments, cc_assignments);

		// Match the timelines

		ArrayList<AliasAssignment> new_timelines = new ArrayList<AliasAssignment>();
		ArrayList<AliasAssignment> revived_timelines = new ArrayList<AliasAssignment>();
		ArrayList<AliasAssignment> stopped_timelines = new ArrayList<AliasAssignment>();
		ArrayList<AliasAssignment> live_timelines = new ArrayList<AliasAssignment>();
		ArrayList<AliasAssignment> dead_timelines = new ArrayList<AliasAssignment>();

		match_timelines_main (db_assignments, cc_assignments,
								new_timelines, revived_timelines, stopped_timelines,
								live_timelines, dead_timelines);

		// Assign timeline IDs to the new timelines

		for (AliasAssignment aa : new_timelines) {
			cc_assignments.set_timeline_id (aa, create_timeline_id (aa.get_primary_id()));
		}

		// Calculate the family time, so it's larger than any of our superseded families

		long family_time = Math.max (sg.task_disp.get_time(), db_assignments.get_family_time() + 1L);

		// Verify invariants

		cc_assignments.check_invariant();
		cc_assignments.check_invariant_db();

		// Begin summary information

		StringBuilder summary_info = new StringBuilder();

		// Report live timelines

		for (AliasAssignment aa : live_timelines) {
			//if (alias_verbose) {
			//	System.out.println ("Alias: Live timeline ID = " + aa.get_timeline_id());
			//}
			summary_info.append (aa.one_line_string() + " live\n");
		}

		// Report dead timelines

		for (AliasAssignment aa : dead_timelines) {
			//if (alias_verbose) {
			//	System.out.println ("Alias: Dead timeline ID = " + aa.get_timeline_id());
			//}
			summary_info.append (aa.one_line_string() + " dead\n");
		}

		// Write tasks for the new timelines

		for (AliasAssignment aa : new_timelines) {

			delete_alias_notification_tasks (aa.get_timeline_id());
		
			OpAliasSplit alias_split_payload = new OpAliasSplit();
			alias_split_payload.setup (family_time);

			PendingTask.submit_task (
				aa.get_timeline_id(),									// event id
				sg.task_sup.get_prompt_exec_time(),						// sched_time
				sg.task_disp.get_time(),								// submit_time
				SUBID_AAFS,												// submit_id
				OPCODE_ALIAS_SPLIT,										// opcode
				0,														// stage
				alias_split_payload.marshal_task());					// details

			//if (alias_verbose) {
			//	System.out.println ("Alias: Split timeline ID = " + aa.get_timeline_id());
			//}
			summary_info.append (aa.one_line_string() + " split\n");
		}

		// Write tasks for the revived timelines

		for (AliasAssignment aa : revived_timelines) {

			delete_alias_notification_tasks (aa.get_timeline_id());
		
			OpAliasRevive alias_revive_payload = new OpAliasRevive();
			alias_revive_payload.setup (family_time);

			PendingTask.submit_task (
				aa.get_timeline_id(),									// event id
				sg.task_sup.get_prompt_exec_time(),						// sched_time
				sg.task_disp.get_time(),								// submit_time
				SUBID_AAFS,												// submit_id
				OPCODE_ALIAS_REVIVE,									// opcode
				0,														// stage
				alias_revive_payload.marshal_task());					// details

			//if (alias_verbose) {
			//	System.out.println ("Alias: Revived timeline ID = " + aa.get_timeline_id());
			//}
			summary_info.append (aa.one_line_string() + " revived\n");
		}

		// Write tasks for the stopped timelines

		for (AliasAssignment aa : stopped_timelines) {

			delete_alias_notification_tasks (aa.get_timeline_id());
		
			OpAliasStop alias_stop_payload = new OpAliasStop();
			alias_stop_payload.setup (family_time);

			PendingTask.submit_task (
				aa.get_timeline_id(),									// event id
				sg.task_sup.get_prompt_exec_time(),						// sched_time
				sg.task_disp.get_time(),								// submit_time
				SUBID_AAFS,												// submit_id
				OPCODE_ALIAS_STOP,										// opcode
				0,														// stage
				alias_stop_payload.marshal_task());						// details

			//if (alias_verbose) {
			//	System.out.println ("Alias: Stopped timeline ID = " + aa.get_timeline_id());
			//}
			summary_info.append (aa.one_line_string() + " stopped\n");
		}

		// Finish the summary information

		summary_info.append ("removed IDs = " + Arrays.toString (cc_assignments.get_all_removed_ids_as_array()));
		String summary_string = summary_info.toString();

		if (alias_verbose) {
			System.out.println ("Alias: Update summary:");
			System.out.println (summary_string);
		}

		sg.log_sup.report_alias_family_updated (family_time, summary_string);

		// Write the new alias family

		AliasFamily.submit_alias_family (null, family_time, cc_assignments);
		return;
	}




	// Get mainshock parameters, given the timeline ID.
	// Parameters:
	//  timeline_id = Timeline ID to search for.
	//  fcmain = Forecast mainshock structure, to be filled in with mainshock parameters.
	// Return values:
	//  RESCODE_SUCCESS - Timeline exists, mainshock parameters are in fcmain and
	//    agree with the current aliases, fcmain.timeline_id contains the timeline ID.
	//  RESCODE_ALIAS_TIMELINE_NOT_FOUND - The timeline ID is not in the alias database.
	//  RESCODE_ALIAS_STOPPED - The timeline ID refers to a stopped timeline.
	// The function throws ComcatException or ComcatConflictException if the operation
	//  cannot be completed due to a problem with Comcat.
	// Any other exception may indicate a problem with the database.
	//
	// Usage notes:
	//
	// This function is used for tasks that contain a timeline ID in the event_id field,
	//  after checking whether a TimelineEntry for that timeline ID exists.
	//
	// If it is a generate task (e.g. forecast), then:
	//  RESCODE_ALIAS_TIMELINE_NOT_FOUND indicates that the alias database is inconsistent
	//    with the timeline database.  The appropriate action is to throw DBCorruptException.
	//  RESCODE_ALIAS_STOPPED indicates that the corresponding event has been deleted
	//    from Comcat.  The appropriate action is to throw ComcatRemovedException, which
	//    invokes the Comcat retry logic to monitor if an event deleted from Comcat is
	//    later restored.
	//
	// Note that a generate task may not call this function if the TimelineEntry does not
	//  exist.  It must fail itself in this case.
	//
	// If it is an intake task (including analyst-intervene and split), then:
	//  RESCODE_ALIAS_TIMELINE_NOT_FOUND indicates that an external user or process
	//    supplied an invalid timeline ID.  The appropriate action is to fail the task.
	//    (It could also indicate a problem with the database, but we ignore this.)
	//  RESCODE_ALIAS_STOPPED indicates that the task refers to an event that is not in
	//    Comcat, or was deleted from Comcat since the original event ID was checked.
	//    The appropriate action is to fail the task.
	//
	// Note that analyst-intervene and intake-sync are not allowed to call this function if
	//  the TimelineEntry exists.  They must operate without calling this function (or any
	//  other function that requires Comcat access).

	public int get_mainshock_for_timeline_id (String timeline_id, ForecastMainshock fcmain) {

		// Make 2 attempts

		for (int attempt = 1; attempt < 3; ++attempt) {

			// Get the alias family for this timeline ID

			AliasFamily alfam = AliasFamily.get_recent_alias_family (0L, 0L, timeline_id, null, null);

			if (alfam == null) {
				return RESCODE_ALIAS_TIMELINE_NOT_FOUND;
			}

			// Read it into an alias assignment list

			AliasAssignmentList db_aalist = alfam.get_assignments();

			// Get the assignment for our timeline

			AliasAssignment db_aa = db_aalist.get_assignment_for_timeline_id (timeline_id);

			if (db_aa == null) {
				throw new DBCorruptException ("AliasSupport.get_mainshock_for_timeline_id: Retrieved alias family does not contain expected timeline ID: " + timeline_id);
			}

			// Assignment list to hold the Comcat results

			AliasAssignmentList cc_aalist = new AliasAssignmentList();

			// If the database assignment has Comcat IDs ...

			if (db_aa.has_comcat_id()) {

				// Query Comcat using the primary ID

				if (!( db_aa.has_primary_id() )) {
					throw new DBCorruptException ("AliasSupport.get_mainshock_for_timeline_id: Retrieved alias family does not have a primary ID: timeline ID: " + timeline_id);
				}

				fcmain.setup_mainshock_poll (db_aa.get_primary_id());

				// If we found the primary ID in Comcat ...

				if (fcmain.mainshock_avail) {

					// Get it into an assignment

					AliasAssignment cc_aa = new AliasAssignment();
					cc_aa.set_comcat_ids_from_array (fcmain.mainshock_id_list);

					// If it matches the database assignment, then we're done

					if (db_aa.is_same_primary_comcat_ids (cc_aa)) {
						fcmain.timeline_id = timeline_id;
						return RESCODE_SUCCESS;
					}

					// Otherwise, add it to the list

					cc_aalist.add_assignment (cc_aa);
				}

				// Otherwise we didn't find the primary ID in Comcat ...

				else {
			
					// Add the primary ID as an absent ID

					cc_aalist.add_absent_id (db_aa.get_primary_id());
				}
			}

			// Otherwise, the database assignment has no Comcat IDs ...

			else {

				// Clear the mainshock parameters available flag

				fcmain.mainshock_avail = false;

				// Loop over removed IDs
				
				for (Iterator<String> it = db_aa.get_removed_id_iterator(); it.hasNext(); ) {
					String removed_id = it.next();

					// If it is not a known Comcat ID ...

					if (!( db_aalist.contains_comcat_id (removed_id) )) {
				
						// Look for it in Comcat

						fcmain.setup_mainshock_poll (removed_id);

						// If it is known in Comcat ...

						if (fcmain.mainshock_avail) {

							// Get it into an assignment

							AliasAssignment cc_aa = new AliasAssignment();
							cc_aa.set_comcat_ids_from_array (fcmain.mainshock_id_list);

							// Add it to the list

							cc_aalist.add_assignment (cc_aa);

							// Exit this loop

							break;
						}

						// Otherwise, add it as an absent ID

						cc_aalist.add_absent_id (removed_id);
					}
				}

				// If we didn't find any of our removed IDs in Comcat, then the timeline is stopped

				if (!( fcmain.mainshock_avail )) {
					return RESCODE_ALIAS_STOPPED;
				}
			}

			// Update the timelines, if this is the first attempt

			if (attempt == 1) {

				if (alias_verbose) {
					System.out.println ("Alias: Update for timeline ID = " + timeline_id);
				}

				update_timelines (db_aalist, cc_aalist);
			}

		}

		// Coming here probably means that Comcat results changed while we were doing this

		throw new ComcatConflictException ("AliasSupport.get_mainshock_for_timeline_id: Comcat data changed during alias update");
	}




	// Get mainshock parameters, given the timeline ID, for a generate task.
	// Parameters:
	//  timeline_id = Timeline ID to search for.
	//  fcmain = Forecast mainshock structure, to be filled in with mainshock parameters.
	// Returns RESCODE_SUCCESS.
	// See documentation above.  This function should only be used by a generate task
	//  (e.g. forecast).  It throws the appropriate exception in case of error, and
	//  so it always returns success.

	public int get_mainshock_for_timeline_id_generate (String timeline_id, ForecastMainshock fcmain) {

		int retval = get_mainshock_for_timeline_id (timeline_id, fcmain);

		if (retval == RESCODE_ALIAS_TIMELINE_NOT_FOUND) {
			throw new DBCorruptException ("AliasSupport.get_mainshock_for_timeline_id_generate: Cannot find timeline ID: " + timeline_id);
		}

		if (retval == RESCODE_ALIAS_STOPPED) {
			throw new ComcatRemovedException ("AliasSupport.get_mainshock_for_timeline_id_generate: Stopped timeline for timeline ID: " + timeline_id);
		}

		return retval;
	}




	// Get mainshock parameters, given the event ID.
	// Parameters:
	//  event_id = Event ID to search for.
	//  fcmain = Forecast mainshock structure, to be filled in with mainshock parameters.
	// Return values:
	//  RESCODE_SUCCESS - Timeline exists, mainshock parameters are in fcmain and
	//    agree with the current aliases, fcmain.timeline_id contains the timeline ID.
	//  RESCODE_ALIAS_EVENT_NOT_IN_COMCAT - The event ID is not known to Comcat.
	//  RESCODE_ALIAS_NEW_EVENT - The event ID is not in the alias database, mainshock
	//    parameters are in fcmain.
	// The function throws ComcatException or ComcatConflictException if the operation
	//  cannot be completed due to a problem with Comcat.
	// Any other exception may indicate a problem with the database.
	//
	// Usage notes:
	//
	// This function is used for tasks that contain a Comcat ID in the event_id field,
	//  and need to convert it into a timeline ID that is used to stage the task.
	//
	// This function can only be used by intake tasks (including analyst-intervene and split).
	//
	// RESCODE_ALIAS_EVENT_NOT_IN_COMCAT indicates that the event ID is not in Comcat,
	//   and no information is available.  (The alias database is not checked.)  The
	//   appropriate action is to fail the task.
	//
	// RESCODE_ALIAS_NEW_EVENT indicates that the event ID is in Comcat, and information
	//   is returned in fcmain, but the event is not in the alias database.  The
	//   appropriate action is to examine the information and determine if it is desired
	//   to continue.  If so, then call write_mainshock_to_new_timeline to insert the
	//   event into the alias database and obtain the timeline ID.

	public int get_mainshock_for_event_id (String event_id, ForecastMainshock fcmain) {

		// Make 2 attempts

		for (int attempt = 1; attempt < 3; ++attempt) {

			// Query Comcat using the event ID

			fcmain.setup_mainshock_poll (event_id);

			// If we didn't find the event in Comcat, return

			if (!( fcmain.mainshock_avail )) {
				return RESCODE_ALIAS_EVENT_NOT_IN_COMCAT;
			}

			// Get it into an assignment

			AliasAssignment cc_aa = new AliasAssignment();
			cc_aa.set_comcat_ids_from_array (fcmain.mainshock_id_list);

			// And an assignment list

			AliasAssignmentList cc_aalist = new AliasAssignmentList();
			cc_aalist.add_assignment (cc_aa);

			// Get the alias family for this set of IDs

			String[] ccids = cc_aa.get_comcat_ids_as_array();

			AliasFamily alfam = AliasFamily.get_recent_alias_family (0L, 0L, null, ccids, null);

			if (alfam == null) {
				return RESCODE_ALIAS_NEW_EVENT;
			}

			// Read it into an alias assignment list

			AliasAssignmentList db_aalist = alfam.get_assignments();

			// Find the assignment for our primary ID (any of our Comcat IDs would do)

			AliasAssignment db_aa = db_aalist.get_assignment_for_comcat_id (cc_aa.get_primary_id());
					
			// If it matches the Comcat assignment, then we're done

			if (db_aa != null) {
				if (db_aa.is_same_primary_comcat_ids (cc_aa)) {
					fcmain.timeline_id = db_aa.get_timeline_id();
					return RESCODE_SUCCESS;
				}
			}

			// Update the timelines, if this is the first attempt

			if (attempt == 1) {

				if (alias_verbose) {
					System.out.println ("Alias: Update for Comcat ID = " + event_id);
				}

				update_timelines (db_aalist, cc_aalist);
			}

		}

		// Coming here probably means that Comcat results changed while we were doing this

		throw new ComcatConflictException ("AliasSupport.get_mainshock_for_event_id: Comcat data changed during alias update");
	}




	// Write mainshock parameters into a new timeline.
	// Parameters:
	//  fcmain = Forecast mainshock structure, containing mainshock parameters.
	// On return, fcmain.timeline_id contains the new timeline name.
	// Any exception may indicate a problem with the database.
	// Note: Before calling this function, you must call get_mainshock_for_event_id and
	// receive a result code of RESCODE_ALIAS_NEW_EVENT.

	public void write_mainshock_to_new_timeline (ForecastMainshock fcmain) {

		// Get the timeline ID to use

		String timeline_id = create_timeline_id (fcmain.mainshock_event_id);

		// Save it to forecast structure for return to caller

		fcmain.timeline_id = timeline_id;

		// Get mainshock into an assignment

		AliasAssignment cc_aa = new AliasAssignment();
		cc_aa.set_timeline_id (timeline_id);
		cc_aa.set_comcat_ids_from_array (fcmain.mainshock_id_list);

		// And an assignment list

		AliasAssignmentList cc_aalist = new AliasAssignmentList();
		cc_aalist.add_assignment (cc_aa);

		// The time

		long family_time = sg.task_disp.get_time();

		// Verify invariants

		cc_aalist.check_invariant();
		cc_aalist.check_invariant_db();

		// Write the new alias family

		sg.log_sup.report_alias_family_created (family_time, cc_aa.one_line_string());

		AliasFamily.submit_alias_family (null, family_time, cc_aalist);
		return;
	}




	// Given a primary ID, find the corresponding timeline ID.
	// Parameters:
	//  primary_id = Primary event ID to search for.
	// Returns the timeline ID, or null if none.
	// Note: This function returns null if the given primary ID is assigned
	// to a timeline but as a secondary ID.
	// Note: This function does not call Comcat nor make any changes to the database,
	// it is strictly a database query.

	public String get_timeline_id_for_primary_id (String primary_id) {

		// Get the alias family for this ID

		String[] ccids = new String[1];
		ccids[0] = primary_id;

		AliasFamily alfam = AliasFamily.get_recent_alias_family (0L, 0L, null, ccids, null);

		if (alfam == null) {
			return null;
		}

		// Read it into an alias assignment list

		AliasAssignmentList db_aalist = alfam.get_assignments();

		// Find the assignment for our primary ID

		AliasAssignment db_aa = db_aalist.get_assignment_for_primary_id (primary_id);

		if (db_aa == null) {
			return null;
		}
	
		// Return the timeline ID

		return db_aa.get_timeline_id();
	}



	
	//----- Construction -----


	// Default constructor.

	public AliasSupport () {}

}

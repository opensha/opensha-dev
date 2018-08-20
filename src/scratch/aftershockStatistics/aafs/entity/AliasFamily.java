package scratch.aftershockStatistics.aafs.entity;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;

import org.bson.types.ObjectId;
import org.mongodb.morphia.annotations.Entity;
import org.mongodb.morphia.annotations.Id;
import org.mongodb.morphia.annotations.Index;
import org.mongodb.morphia.annotations.Indexed;
import org.mongodb.morphia.annotations.IndexOptions;
import org.mongodb.morphia.annotations.Indexes;
import org.mongodb.morphia.annotations.Field;
import org.mongodb.morphia.annotations.Transient;
import org.mongodb.morphia.annotations.Embedded;

import org.mongodb.morphia.Datastore;
import org.mongodb.morphia.FindAndModifyOptions;

import org.mongodb.morphia.query.Query;
import org.mongodb.morphia.query.UpdateOperations;
import org.mongodb.morphia.query.MorphiaIterator;

import scratch.aftershockStatistics.aafs.AliasAssignment;
import scratch.aftershockStatistics.aafs.AliasAssignmentList;
import scratch.aftershockStatistics.aafs.MongoDBUtil;
import scratch.aftershockStatistics.aafs.RecordKey;
import scratch.aftershockStatistics.aafs.RecordPayload;
import scratch.aftershockStatistics.aafs.RecordIterator;

import scratch.aftershockStatistics.util.MarshalImpArray;
import scratch.aftershockStatistics.util.MarshalImpJsonReader;
import scratch.aftershockStatistics.util.MarshalImpJsonWriter;
import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;




/**
 * Holds an AAFS alias family in the MongoDB database.
 * Author: Michael Barall 06/23/2018.
 *
 * The collection "alias" holds the alias families.
 *
 * The purpose of this class is to contain the mapping from Comcat IDs to timeline IDs.
 * This is relatively complicated because each timeline ID can correspond to several Comcat IDs;
 * the set of Comcat IDs associated with a timeline ID can change over time;  the mapping has
 * to be inferred from id lists returned by Comcat;  events can be split (Comcat IDs that refer
 * to the same event can later refer to different events);  and events can be merged (Comcat IDs
 * that refer to different events can later refer to the same event).
 *
 * Introduce some terminology.  Call two Comcat IDs "siblings" if at some time they were
 * associated to the same event.  (A Comcat ID is considered to be its own sibling.)
 * Call two Comcat IDs C1 and CN "cousins" if there are Comcat IDs C2, C3, ... such that
 * C1 and C2 are siblings, and C2 and C3 are siblings, etc.  Then the "cousin" relation is
 * the transitive closure of the "sibling" relation and hence is an equivalence relation.
 * An equivalence class of the "cousin" relation is called a "family".
 *
 * Each object of this class contains:
 * -- All the Comcat IDs making up a family at a given time.
 * -- All the timeline IDs that have ever been associated with any of the Comcat IDs.
 * -- A time stamp.
 * -- A map, which specifies which of the timeline IDs is currently associated to each Comcat ID.
 * -- For each timeline ID with associated Comcat IDs, an indication of which Comcat ID is primary.
 * (It is possible for a Comcat ID to not be associated with any timeline ID, e.g., if it has
 * been deleted from Comcat.  It is also possible for a timeline ID to not be associated with
 * any Comcat ID, e.g., as a result of merging.)
 *
 * If object O1 has a later time stamp than object O2, and if they have at least one Comcat ID
 * in common, then we say that O1 "supersedes" O2.  If O1 supersedes O2 then it follows that
 * O1 contains all the Comcat IDs that are in O2, that is, O1 is a superset of O2.  We say that
 * an object O is "active" if it is not superseded by any other object.  Because of the superset
 * property, a given Comcat ID will appear in exactly one active object, which can be found by
 * searching for the most recent object that contains the Comcat ID.
 *
 * By definition, the Comcat IDs assocated with a given timeline ID at a given time are all
 * siblings (and hence cousins).  We impose the following continuity condition: For a given
 * timeline ID T, whenever the set of Comcat IDs associated with T changes, the new set must
 * either be empty or else contain at least one Comcat ID that was previously associated with T.
 * The continuity condition implies that all the Comcat IDs ever associated with T are cousins,
 * and hence lie in a single family.  So a given timeline ID will appear in exactly one active
 * object, which contains all the Comcat IDs ever associated to that object, and can be found
 * by searching for the most recent object that contains the timeline ID.
 *
 * This structure is designed so that changes to the mapping can be effected with a single
 * write to the database, eliminating the need for any sort of transactions.  Also, the mapping
 * in effect at any time in the past can be recovered by disregarding any more recent objects.
 *
 * It should be noted that the most common case is for the set of Comcat IDs associated with
 * a given event to never change.  In this common case, the corresponding object contains a
 * single timeline ID and is never superseded.
 */
@Entity(value = "alias", noClassnameStored = true)
@Indexes({
	@Index(fields = {@Field("timeline_ids")}, options = @IndexOptions(name = "famtlid")),
	@Index(fields = {@Field("comcat_ids")}, options = @IndexOptions(name = "famccid")),
	@Index(fields = {@Field("family_time")}, options = @IndexOptions(name = "famtime"))
})
public class AliasFamily implements java.io.Serializable {

	//----- Envelope information -----

	// Globally unique identifier for this alias family.
	// This is the MongoDB identifier.
	// Note that ObjectId implements java.io.Serializable.
	// This is set to the same value as the id of the task that generated the alias family.

    @Id
    private ObjectId id;

	// Time that this alias family was created, in milliseconds since the epoch.
	// The collection is indexed on this field, so alias families in a given
	// time span can be obtained with a database query efficiently.
	// In practice, every family in the collection must have a different value of family_time
	// (so time values may need to be altered slightly), strictly monotonically increasing
	// in the order that entries are written.

	//@Indexed(options = @IndexOptions(name = "famtime"))
	private long family_time;

	//----- Action information -----

	// List of timeline IDs for this family.
	// This cannot be null or empty, and the array elements cannot be null.
	// This can be used to find all alias families pertaining to a timeline.

	//@Indexed(options = @IndexOptions(name = "famtlid"))
	private String[] timeline_ids;

	// List of Comcat IDs for this family.
	// This cannot be null or empty, and the array elements cannot be null.
	// This can be used to find all alias families pertaining to one or more Comcat IDs.
	// The ordering of Comcat IDs is significant.  First come the current Comcat IDs for the
	// first timeline, with the primary ID first.  Then come the current Comcat IDs for the
	// second timeline, with the primary ID first.  And so on.  After all the current Comcat
	// IDs for all the timelines, then come the Comcat IDs that were previously associated
	// with one or more timelines but are not currently associated with any timeline.

	//@Indexed(options = @IndexOptions(name = "famccid"))
	private String[] comcat_ids;

	// Encoded list of bindings for this alias family.

	private int[] enc_bindings;




	//----- Getters and setters -----

    private ObjectId get_id() {
        return id;
    }

    private void set_id (ObjectId id) {
        this.id = id;
    }

    public long get_family_time() {
        return family_time;
    }

    private void set_family_time (long family_time) {
        this.family_time = family_time;
    }

    public String[] get_timeline_ids() {
        return timeline_ids.clone();
    }

    private void set_timeline_ids (String[] timeline_ids) {
        this.timeline_ids = timeline_ids.clone();
    }

    public String[] get_comcat_ids() {
        return comcat_ids.clone();
    }

    private void set_comcat_ids (String[] comcat_ids) {
        this.comcat_ids = comcat_ids.clone();
    }

    public int[] get_enc_bindings() {
        return enc_bindings.clone();
    }

    private void set_enc_bindings (int[] enc_bindings) {
        this.enc_bindings = enc_bindings.clone();
    }




	// toString - Convert to string.

	@Override
	public String toString() {
		String str = "AliasFamily\n"
			+ "\tid: " + ((id == null) ? ("null") : (id.toHexString())) + "\n"
			+ "\tfamily_time: " + family_time + "\n"
			+ "\ttimeline_ids: " + Arrays.toString (timeline_ids) + "\n"
			+ "\tcomcat_ids: " + Arrays.toString (comcat_ids) + "\n"
			+ "\tenc_bindings: " + Arrays.toString (enc_bindings);
		return str;
	}




	/**
	 * get_record_key - Get the record key for this alias family.
	 */
	public RecordKey get_record_key () {
		return new RecordKey(id);
	}




	/**
	 * set_record_key - Set the record key for this alias family.
	 * @param key = Record key. Can be null.
	 */
	private void set_record_key (RecordKey key) {
		if (key == null) {
			id = null;
		} else {
			id = key.getId();
		}
		return;
	}




	// Set the timeline IDs, Comcat IDs, and bindings from a list of assignments.
	// The family time in the assigment list is not used.

	private void set_assignments (AliasAssignmentList assignments) {

		if (!( assignments != null )) {
			throw new IllegalArgumentException("AliasFamily.set_assigments: No assignment list supplied");
		}

		// Get the required lengths of our arrays

		int timeline_ids_length = assignments.get_assignment_count();
		if (!( timeline_ids_length > 0 )) {
			throw new IllegalArgumentException("AliasFamily.set_assigments: No timeline IDs supplied");
		}

		int comcat_ids_length = assignments.get_all_comcat_ids_count() + assignments.get_all_removed_ids_count();
		if (!( comcat_ids_length > 0 )) {
			throw new IllegalArgumentException("AliasFamily.set_assigments: No Comcat IDs supplied");
		}

		int enc_bindings_length = 2*timeline_ids_length + assignments.get_total_removed_ids_count();

		// Allocate the arrays

		timeline_ids = new String[timeline_ids_length];
		comcat_ids = new String[comcat_ids_length];
		enc_bindings = new int[enc_bindings_length];

		// Indexes into the arrays

		int tlix = 0;
		int ccix = 0;
		int ebix = 2*timeline_ids_length;

		// Maps that invert the timeline and Comcat arrays

		HashMap<String, Integer> tlmap = new HashMap<String, Integer>();
		HashMap<String, Integer> ccmap = new HashMap<String, Integer>();

		// This array will hold the list of removed IDs for each timeline, in the same slots as enc_bindings

		String[] removed_ids = new String[enc_bindings_length];

		// Loop over assignments

		for (Iterator<AliasAssignment> it = assignments.get_assigment_iterator(); it.hasNext(); ) {
			AliasAssignment assignment = it.next();

			// Get the timeline ID

			timeline_ids[tlix] = assignment.get_timeline_id();
			if (!( timeline_ids[tlix] != null )) {
				throw new IllegalArgumentException("AliasFamily.set_assigments: Found assignment with no timeline ID");
			}

			// Insert into reverse map and check for duplicate

			if (tlmap.put (timeline_ids[tlix], new Integer(tlix)) != null) {
				throw new IllegalArgumentException("AliasFamily.set_assigments: Found duplicate timeline ID : " + timeline_ids[tlix]);
			}

			// Get the list of Comcat IDs

			ccix = assignment.get_comcat_ids_as_array (comcat_ids, ccix);

			// Get the list of removed IDs

			ebix = assignment.get_removed_ids_as_array (removed_ids, ebix);

			// Save end-of-list indexes

			enc_bindings[tlix] = ccix;
			enc_bindings[tlix + timeline_ids_length] = ebix;

			// Next timeline

			++tlix;
		}

		// Now append the list of all removed IDs

		ccix = assignments.get_all_removed_ids_as_array (comcat_ids, ccix);

		// Check final Indexes

		if (!( tlix == timeline_ids_length && ccix == comcat_ids_length && ebix == enc_bindings_length )) {
			throw new IllegalArgumentException("AliasFamily.set_assigments: Element count mismatch");
		}

		// Construct the reverse map of Comcat IDs and check for duplicates

		for (ccix = 0; ccix < comcat_ids_length; ++ccix) {
			if (ccmap.put (comcat_ids[ccix], new Integer(ccix)) != null) {
				throw new IllegalArgumentException("AliasFamily.set_assigments: Found duplicate Comcat ID : " + comcat_ids[ccix]);
			}
		}

		// Convert the removed IDs to indexes

		for (ebix = 2*timeline_ids_length; ebix < enc_bindings_length; ++ebix) {
			enc_bindings[ebix] = ccmap.get (removed_ids[ebix]);
		}

		return;
	}




	// Get the timeline IDs, Comcat IDs, and bindings as a list of assignments.
	// Also sets the family time in the assigment list.

	public AliasAssignmentList get_assignments () {

		AliasAssignmentList assignments = new AliasAssignmentList();

		// Get the number of timelines

		int timeline_ids_length = timeline_ids.length;

		// Indexes into the arrays

		int ccix = 0;
		int ebix = 2*timeline_ids_length;

		// Loop over timelines

		for (int tlix = 0; tlix < timeline_ids_length; ++tlix) {

			AliasAssignment assignment = new AliasAssignment();

			// Set the timeline ID

			assignment.set_timeline_id (timeline_ids[tlix]);

			// Set the Comcat IDs

			int hi = enc_bindings[tlix];
			assignment.set_comcat_ids_from_array (comcat_ids, ccix, hi);
			ccix = hi;

			// Set the removed IDs

			hi = enc_bindings[tlix + timeline_ids_length];
			for ( ; ebix < hi; ++ebix) {
				assignment.add_removed_id (comcat_ids[enc_bindings[ebix]]);
			}

			// Add to the list

			assignments.add_assignment (assignment);
		}

		// Set the time stamp

		assignments.set_family_time (family_time);

		return assignments;
	}




	/**
	 * submit_alias_family - Submit an alias family.
	 * @param key = Record key associated with this alias family. Can be null to assign a new one.
	 * @param family_time = Time of this alias family, in milliseconds
	 *                      since the epoch. Must be positive.
	 * @param timeline_ids = List of timeline IDs associated with this alias family. Cannot be null or empty.
	 * @param comcat_ids = List of Comcat IDs associated with this alias family. Cannot be null or empty.
	 * @param enc_bindings = Encoded bindings for the list of assignments. Cannot be null or empty.
	 * @return
	 * Returns the new family.
	 */
	private static AliasFamily submit_alias_family (RecordKey key, long family_time, String[] timeline_ids,
			String[] comcat_ids, int[] enc_bindings) {

		// Check conditions

		if (!( family_time > 0L
			&& timeline_ids != null
			&& timeline_ids.length > 0
			&& comcat_ids != null
			&& comcat_ids.length > 0
			&& enc_bindings != null
			&& enc_bindings.length >= 2*timeline_ids.length )) {
			throw new IllegalArgumentException("AliasFamily.submit_alias_family: Invalid alias family parameters");
		}

		for (String timeline_id : timeline_ids) {
			if (!( timeline_id != null )) {
				throw new IllegalArgumentException("AliasFamily.submit_alias_family: Invalid alias family parameters");
			}
		}

		for (String comcat_id : comcat_ids) {
			if (!( comcat_id != null )) {
				throw new IllegalArgumentException("AliasFamily.submit_alias_family: Invalid alias family parameters");
			}
		}

		// Construct the alias family object

		AliasFamily alfam = new AliasFamily();
		alfam.set_record_key (key);
		alfam.set_family_time (family_time);
		alfam.set_timeline_ids (timeline_ids);
		alfam.set_comcat_ids (comcat_ids);
		alfam.set_enc_bindings (enc_bindings);

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(alfam);
		
		return alfam;
	}




	/**
	 * submit_alias_family - Submit an alias family.
	 * @param key = Record key associated with this alias family. Can be null to assign a new one.
	 * @param family_time = Time of this alias family, in milliseconds
	 *                      since the epoch. Must be positive.
	 * @param assignments = List of assignments for this alias family. Cannot be null or empty.
	 * @return
	 * Returns the new family.
	 */
	public static AliasFamily submit_alias_family (RecordKey key, long family_time, AliasAssignmentList assignments) {

		// Check conditions

		if (!( family_time > 0L
			&& assignments != null )) {
			throw new IllegalArgumentException("AliasFamily.submit_alias_family: Invalid alias family parameters");
		}

		// Construct the alias family object

		AliasFamily alfam = new AliasFamily();
		alfam.set_record_key (key);
		alfam.set_family_time (family_time);
		alfam.set_assignments (assignments);

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(alfam);
		
		return alfam;
	}




	/**
	 * store_alias_family - Store an alias family into the database.
	 * This is primarily for restoring from backup.
	 */
	public static AliasFamily store_alias_family (AliasFamily alfam) {

		// Call MongoDB to store into database

		Datastore datastore = MongoDBUtil.getDatastore();
		datastore.save(alfam);
		
		return alfam;
	}




	/**
	 * get_alias_family_for_key - Get the alias family with the given key.
	 * @param key = Record key. Cannot be null or empty.
	 * Returns the alias family, or null if not found.
	 */
	public static AliasFamily get_alias_family_for_key (RecordKey key) {

		if (!( key != null && key.getId() != null )) {
			throw new IllegalArgumentException("AliasFamily.get_alias_family_for_key: Missing or empty record key");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query: Select id == key

		Query<AliasFamily> query = datastore.createQuery(AliasFamily.class)
											.filter("id ==", key.getId());

		// Run the query

		AliasFamily alfam = query.get();

		return alfam;
	}




	/**
	 * get_alias_family_range - Get a range of alias families, reverse-sorted by action time.
	 * @param family_time_lo = Minimum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no minimum.
	 * @param family_time_hi = Maximum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no maximum.
	 * @param timeline_id = Timeline id. Can be null to return entries for all timelines.
	 * @param comcat_ids = Comcat id list. Can be null or empty to return entries for all Comcat ids.
	 *                     If specified, return entries associated with any of the given ids.
	 * @param family_time_div_rem = 2-element array containing divisor (element 0) and remainder (element 1) for
	 *                              action time modulus. Can be null, or contain zeros, for no modulus test.
	 */
	public static List<AliasFamily> get_alias_family_range (long family_time_lo, long family_time_hi, String timeline_id, String[] comcat_ids, long[] family_time_div_rem) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<AliasFamily> query = datastore.createQuery(AliasFamily.class);

		// Select by timeline_id

		if (timeline_id != null) {
			query = query.filter("timeline_ids ==", timeline_id);
		}

		// Select by comcat_ids

		if (comcat_ids != null) {
			if (comcat_ids.length > 0) {
				query = query.filter("comcat_ids in", comcat_ids);
			}
		}

		// Select entries with family_time >= family_time_lo

		if (family_time_lo > 0L) {
			query = query.filter("family_time >=", new Long(family_time_lo));
		}

		// Select entries with family_time <= family_time_hi

		if (family_time_hi > 0L) {
			query = query.filter("family_time <=", new Long(family_time_hi));
		}

		// Select entries with family_time % family_time_div_rem[0] == family_time_div_rem[1]

		if (family_time_div_rem != null) {
			if (family_time_div_rem[0] > 0L) {
				Long[] div_rem = new Long[2];
				div_rem[0] = new Long(family_time_div_rem[0]);
				div_rem[1] = new Long(family_time_div_rem[1]);
				query = query.filter("family_time mod", div_rem);
			}
		}

		// Sort by family_time in descending order (most recent first)

		query = query.order("-family_time");

		// Run the query

		List<AliasFamily> entries = query.asList();

		return entries;
	}




	/**
	 * fetch_alias_family_range - Iterate a range of alias families, reverse-sorted by action time.
	 * @param family_time_lo = Minimum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no minimum.
	 * @param family_time_hi = Maximum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no maximum.
	 * @param timeline_id = Timeline id. Can be null to return entries for all timelines.
	 * @param comcat_ids = Comcat id list. Can be null or empty to return entries for all Comcat ids.
	 *                     If specified, return entries associated with any of the given ids.
	 * @param family_time_div_rem = 2-element array containing divisor (element 0) and remainder (element 1) for
	 *                              action time modulus. Can be null, or contain zeros, for no modulus test.
	 */
	public static RecordIterator<AliasFamily> fetch_alias_family_range (long family_time_lo, long family_time_hi, String timeline_id, String[] comcat_ids, long[] family_time_div_rem) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<AliasFamily> query = datastore.createQuery(AliasFamily.class);

		// Select by timeline_id

		if (timeline_id != null) {
			query = query.filter("timeline_ids ==", timeline_id);
		}

		// Select by comcat_ids

		if (comcat_ids != null) {
			if (comcat_ids.length > 0) {
				query = query.filter("comcat_ids in", comcat_ids);
			}
		}

		// Select entries with family_time >= family_time_lo

		if (family_time_lo > 0L) {
			query = query.filter("family_time >=", new Long(family_time_lo));
		}

		// Select entries with family_time <= family_time_hi

		if (family_time_hi > 0L) {
			query = query.filter("family_time <=", new Long(family_time_hi));
		}

		// Select entries with family_time % family_time_div_rem[0] == family_time_div_rem[1]

		if (family_time_div_rem != null) {
			if (family_time_div_rem[0] > 0L) {
				Long[] div_rem = new Long[2];
				div_rem[0] = new Long(family_time_div_rem[0]);
				div_rem[1] = new Long(family_time_div_rem[1]);
				query = query.filter("family_time mod", div_rem);
			}
		}

		// Sort by family_time in descending order (most recent first)

		query = query.order("-family_time");

		// Run the query

		MorphiaIterator<AliasFamily, AliasFamily> morphia_iterator = query.fetch();

		return new RecordIterator<AliasFamily>(morphia_iterator);
	}




	/**
	 * get_recent_alias_family - Get the most recent in a range of alias families.
	 * @param family_time_lo = Minimum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no minimum.
	 * @param family_time_hi = Maximum action time, in milliseconds since the epoch.
	 *                         Can be 0L for no maximum.
	 * @param timeline_id = Timeline id. Can be null to return entries for all timelines.
	 * @param comcat_ids = Comcat id list. Can be null or empty to return entries for all Comcat ids.
	 *                     If specified, return entries associated with any of the given ids.
	 * @param family_time_div_rem = 2-element array containing divisor (element 0) and remainder (element 1) for
	 *                              action time modulus. Can be null, or contain zeros, for no modulus test.
	 * Returns the matching alias family with the greatest family_time (most recent),
	 * or null if there is no matching alias family.
	 */
	public static AliasFamily get_recent_alias_family (long family_time_lo, long family_time_hi, String timeline_id, String[] comcat_ids, long[] family_time_div_rem) {

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Construct the query

		Query<AliasFamily> query = datastore.createQuery(AliasFamily.class);

		// Select by timeline_id

		if (timeline_id != null) {
			query = query.filter("timeline_ids ==", timeline_id);
		}

		// Select by comcat_ids

		if (comcat_ids != null) {
			if (comcat_ids.length > 0) {
				query = query.filter("comcat_ids in", comcat_ids);
			}
		}

		// Select entries with family_time >= family_time_lo

		if (family_time_lo > 0L) {
			query = query.filter("family_time >=", new Long(family_time_lo));
		}

		// Select entries with family_time <= family_time_hi

		if (family_time_hi > 0L) {
			query = query.filter("family_time <=", new Long(family_time_hi));
		}

		// Select entries with family_time % family_time_div_rem[0] == family_time_div_rem[1]

		if (family_time_div_rem != null) {
			if (family_time_div_rem[0] > 0L) {
				Long[] div_rem = new Long[2];
				div_rem[0] = new Long(family_time_div_rem[0]);
				div_rem[1] = new Long(family_time_div_rem[1]);
				query = query.filter("family_time mod", div_rem);
			}
		}

		// Sort by family_time in descending order (most recent first)

		query = query.order("-family_time");

		// Run the query

		AliasFamily alfam = query.get();

		return alfam;
	}




	/**
	 * delete_alias_family - Delete an alias family.
	 * @param alfam = Existing alias family to delete.
	 * @return
	 */
	public static void delete_alias_family (AliasFamily alfam) {

		// Check conditions

		if (!( alfam != null && alfam.get_id() != null )) {
			throw new IllegalArgumentException("AliasFamily.delete_alias_family: Invalid parameters");
		}

		// Get the MongoDB data store

		Datastore datastore = MongoDBUtil.getDatastore();

		// Run the delete

		datastore.delete(alfam);
		
		return;
	}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 11001;

	private static final String M_VERSION_NAME = "AliasFamily";

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

		String sid = id.toHexString();
		writer.marshalString      ("id"          , sid         );
		writer.marshalLong        ("family_time" , family_time );
		writer.marshalStringArray ("timeline_ids", timeline_ids);
		writer.marshalStringArray ("comcat_ids"  , comcat_ids  );
		writer.marshalIntArray    ("enc_bindings", enc_bindings);
	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents

		String sid;
		sid          = reader.unmarshalString      ("id"          );
		family_time  = reader.unmarshalLong        ("family_time" );
		timeline_ids = reader.unmarshalStringArray ("timeline_ids");
		comcat_ids   = reader.unmarshalStringArray ("comcat_ids"  );
		enc_bindings = reader.unmarshalIntArray    ("enc_bindings");
		id = new ObjectId(sid);

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

	public AliasFamily unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}


}

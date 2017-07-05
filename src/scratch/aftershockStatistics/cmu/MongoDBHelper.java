package scratch.aftershockStatistics.cmu;

import java.net.UnknownHostException;
import java.util.Arrays;

import org.mongodb.morphia.Datastore;
import org.mongodb.morphia.Morphia;

import com.mongodb.*;

public enum MongoDBHelper {
    INSTANCE;

    private DB db;
    private Datastore datastore;


    private MongoDBHelper() {
        Configuration config = new Configuration();

        try {

            MongoCredential credentials = MongoCredential.createCredential(config.getDb_user(), config.getDb_name(), config.getDb_password().toCharArray());
            ServerAddress serverAddress = new ServerAddress(config.getDb_host(), config.getDb_port());
            MongoClient mongoClient = new MongoClient(serverAddress , Arrays.asList(credentials));

            this.db = mongoClient.getDB(config.getDb_name());

            Morphia morphia = new Morphia();

            this.datastore = morphia.createDatastore(mongoClient, config.getDb_name());

            morphia.mapPackage("package");
        } catch (UnknownHostException e) {
            e.printStackTrace();
        }

    }

    public DB getDB() {
        return this.db;
    }

    public Datastore getDatastore() {
        return this.datastore;
    }
}
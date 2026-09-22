class Utils {
    static String sanitize(String name) {
        return name.replaceAll('[^A-Za-z0-9_-]', '_').toLowerCase()
    }

    //
    // Whether a meta value is set. null, blank strings and empty lists
    // (nf-schema's value for an empty samplesheet column) all count as unset.
    //
    static boolean hasValue(value) {
        if( value instanceof Collection )
            return value.any { hasValue(it) }
        return value != null && value.toString().trim()
    }
}
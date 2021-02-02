<template>
  <div id="wrapper">
    <div class="table-responsive">
      <div v-if="showWarning" class="alert alert-warning">
        <strong
          >This function does not correspond to the selected 
          open reading frame.</strong
        >
      </div>
      <table id="blast-table" class="table table-hover">
        <thead>
          <tr>
            <th scope="col">Accession</th>
            <th width="300px" scope="col">Description</th>
            <th scope="col">Evalue</th>
            <th scope="col">% Identity</th>
            <th scope="col">Query</th>
            <th scope="col">Subject</th>
            <th scope="col">Action</th>
          </tr>
        </thead>
        <tbody>
          <tr v-for="alignment in blastResults" :key="alignment.accession">
            <td>{{ alignment.accession }}</td>
            <td>{{ alignment.title }}</td>
            <td>{{ alignment.evalue }}</td>
            <td>{{ alignment.percent_identity }}%</td>
            <td>{{ alignment.query_from }} - {{ alignment.query_to }}</td>
            <td>{{ alignment.hit_from }} - {{ alignment.hit_to }}</td>
            <td>
              <button
                v-if="allowSelect"
                class="btn btn-dark btn-sm"
                @click="setFunction(alignment.title, alignment.accession)"
              >
                <strong>Select</strong>
              </button>
              <button
                v-if="!allowSelect"
                class="btn btn-dark btn-sm"
                @click="showWarning = true"
              >
                <strong>Select</strong>
              </button>
            </td>
          </tr>
        </tbody>
      </table>
      <div>
        <h4 v-if="blastResults.length === 0">
          <em>No hits found.</em>
        </h4>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: "BlastResults",
  props: {
    blastResults: Array,
    allowSelect: false,
    showWarning: false,
  },

  data() {
    return {};
  },

  methods: {

    /**
     * Sets the CDS function from the blast results.
     * @param {string} funct the function to be added.
     */
    setFunction(funct, accession) {
      let result = funct.match(/(.*) \[.*\]/);
      this.$emit("newFunction", result[1] + '##' + accession);
    },

  },
};
</script>

<style scoped>
/* ----- Blast Table ----- */
.table-responsive {
  max-height: 500px;
  overflow-y: auto;
  display: inline-block;
}

#blast-table {
  text-align: left;
}

.table-responsive thead th {
  position: sticky;
  top: 0;
  background: #eee;
  border: darkgray;
}
</style>
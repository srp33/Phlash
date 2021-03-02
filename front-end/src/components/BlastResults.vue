<template>
  <div id="wrapper">
    <div class="table-responsive">
      <div v-if="showWarning" class="alert alert-warning">
        <strong
          >This function does not correspond to the selected open reading
          frame.</strong
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
            <td style="font-size: 1.4em">{{ alignment.accession }}</td>
            <td>{{ alignment.title }}</td>
            <td style="font-size: 1.4em">{{ alignment.evalue }}</td>
            <td style="font-size: 1.4em">{{ alignment.percent_identity }}%</td>
            <td style="font-size: 1.4em">
              {{ alignment.query_from }} - {{ alignment.query_to }}
            </td>
            <td style="font-size: 1.4em">
              {{ alignment.hit_from }} - {{ alignment.hit_to }}
            </td>
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
        <h4 v-if="this.blastResults.length === 0">
          <em>No hits found.</em>
        </h4>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'BlastResults',
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
      this.$emit('newFunction', result[1] + '##' + accession);
    },
  },
};
</script>

<style scoped>
.table-responsive {
  max-height: 35em;
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
  background-color: white;
  border: darkgray;
  font-size: 1.4em;
}

.btn-dark {
  font-size: 15pt;
}
</style>
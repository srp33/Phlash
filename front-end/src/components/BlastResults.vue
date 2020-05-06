<template>
  <div id="wrapper">
    <div class="table-responsive">
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
              <button class="btn btn-dark btn-sm" @click="setFunction(alignment.title)">
                <strong>Select</strong>
              </button>
            </td>
          </tr>
        </tbody>
      </table>
      <div v-if="blastResults.length===0">
        <h4>
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
    blastResults: Array
  },
  data() {
    return {
    };
  },
  methods: {
    setFunction(funct) {
      let result = funct.match(/(.*)\[.*\]/);
      this.$emit("newFunction", result[1]);
    }
  }
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
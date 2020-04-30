<template>
  <div class="container">
    <loading loader="dots" :active.sync="isLoading" :is-full-page="true" :height="200" :width="200" color="#2279b6"></loading>
    <div class="headers">
      <h1>ID: {{ $route.params.cdsID }}</h1>
      <h4>Start: {{ currentCDS.start }}</h4>
      <h4>Stop: {{ currentCDS.stop }}</h4>
      <h4>Status: {{ currentCDS.status }}</h4>
    </div>
    <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>
        Choose a new start for this gene call based on the information given
        below or keep the current start.
      </p>
      <p><strong>Your selected start position:</strong> {{ newStart }}</p>
      <p><strong>Your selected function:</strong> {{ newFunction }}</p>
      <button type="button" class="btn btn-light" @click="editCDS" v-if="newFunction !== 'None' && newStart !== 0">
        <strong>Update {{ $route.params.cdsID }}</strong>
      </button>
    </div>
    <div class="columns-wrapper">
      <h3 style="text-align: center; margin: 40px;">GeneMark's Coding Potential</h3>
    </div>
    <div class="coding-potential-table">
      <h5 style="text-align: center; margin: 20px;">Average Coding Potential Per Frame</h5>
      <div class="table-responsive">
        <table id="cp-table" class="table table-hover">
          <caption>Frames are counted 1-6 (direct 1-3 and complementary 4-6).</caption>
          <thead>
            <tr>
              <th scope="col">Start Position</th>
              <th scope="col">Frame 1</th>
              <th scope="col">Frame 2</th>
              <th scope="col">Frame 3</th>
              <th scope="col">Frame 4</th>
              <th scope="col">Frame 5</th>
              <th scope="col">Frame 6</th>
              <th scope="col">Action</th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(data, start, index) in startOptions" :key="index">
              <th v-if="start===currentCDS.start.toString()">{{ start }} (current)</th>
              <th v-else>{{ start }}</th>
              <td>{{ data.frame_1 }}</td>
              <td>{{ data.frame_2 }}</td>
              <td>{{ data.frame_3 }}</td>
              <td>{{ data.frame_4 }}</td>
              <td>{{ data.frame_5 }}</td>
              <td>{{ data.frame_6 }}</td>
              <td>
                <button class="btn btn-dark btn-sm" @click="setStart(start)">
                  <strong>Select</strong>
                </button>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
    <div class="coding-potential-graphs">
      <h5 style="text-align: center; margin: 20px;">Coding Potential Per Frame</h5>
      <div v-if="dataExists">
        <Graphs :data1="data1" :data2="data2" :data3="data3" :data4="data4" :data5="data5" :data6="data6" :start="currentCDS.start" :stop="currentCDS.stop" />
      </div>
      <div v-else>
        <p>Loading graphs...</p>
      </div>
    </div>
    <hr />
    <div class="blast-results">
      <h3 style="text-align: center; margin: 20px;">BLAST Results</h3>
      <div id="accordion">
        <div class="card" v-for="(item, key) in blastResults" :key="key">
          <div class="card-header">
            <h4 class="mb-0">
              <button class="btn btn-light btn-blast" data-toggle="collapse" aria-expanded="false" v-bind:data-target="'#'+key" v-bind:aria-controls="key">
                <strong>{{ key }}</strong>
              </button>
            </h4>
          </div>
          <div v-bind:id="key" class="collapse" data-parent="#accordion">
            <div class="card-body">
              <BlastResults :blastResults="item" @newFunction="setFunction" />
            </div>
          </div>
        </div>
      </div>
    </div>
    <div class="info-bottom">
      <p><strong>Your selected start position:</strong> {{ newStart }}</p>
      <p><strong>Your function selection:</strong> {{ newFunction }}</p>
      <button type="button" class="btn btn-light" @click="editCDS" v-if="newFunction !== 'None' && newStart !== 0">
        <strong>Update {{ $route.params.cdsID }}</strong>
      </button>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import BlastResults from "../components/BlastResults.vue";
import Graphs from "../components/Graphs.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Fail",
  components: {
    BlastResults,
    Graphs,
    Loading
  },
  data() {
    return {
      startOptions: [],
      currentCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: ""
      },
      updatedCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: ""
      },
      blastResults: [],
      showFunction: false,
      showStart: false,
      newFunction: "None",
      newStart: 0,
      dataExists: false,
      isLoading: true,
      data1: [{ x: [], y: [] }],
      data2: [{ x: [], y: [] }],
      data3: [{ x: [], y: [] }],
      data4: [{ x: [], y: [] }],
      data5: [{ x: [], y: [] }],
      data6: [{ x: [], y: [] }]
    };
  },
  created() {
    this.getData(this.$route.params.cdsID);
  },
  methods: {
    getData(cdsID) {
      axios.get(`http://localhost:5000/api/annotations/fail/${this.$route.params.phageID}/${cdsID}`)
        .then(response => {
          this.currentCDS = response.data.cds;
          this.blastResults = response.data.blast;
          this.startOptions = response.data.start_options;
          this.data1 = [{
            x: response.data.x_data,
            y: response.data.y_data_1
          }];
          this.data2 = [{
            x: response.data.x_data,
            y: response.data.y_data_2
          }];
          this.data3 = [{
            x: response.data.x_data,
            y: response.data.y_data_3
          }];
          this.data4 = [{
            x: response.data.x_data,
            y: response.data.y_data_4
          }];
          this.data5 = [{
            x: response.data.x_data,
            y: response.data.y_data_5
          }];
          this.data6 = [{
            x: response.data.x_data,
            y: response.data.y_data_6
          }];
          this.dataExists = true;
          this.isLoading = false;
        })
        .catch(error => {
          console.error(error);
        });
    },
    editCDS() {
      this.updatedCDS = this.currentCDS;
      this.updatedCDS.start = this.newStart;
      this.updatedCDS.function = this.newFunction;
      const payload = {
        id: this.updatedCDS.id,
        start: this.updatedCDS.start,
        stop: this.updatedCDS.stop,
        strand: this.updatedCDS.strand,
        function: this.updatedCDS.function,
        status: "Pass"
      };
      this.updateCDS(payload, this.updatedCDS.id);
    },
    updateCDS(payload, cdsID) {
      axios.put(`http://localhost:5000/api/annotations/fail/${this.$route.params.phageID}/${cdsID}`,
          payload
        )
        .then(() => {
          this.$router.push(`/annotations/${this.$route.params.phageID}`);
        })
        .catch(error => {
          console.error(error);
        });
    },
    keepOriginal() {
      const payload = {
        id: this.currentCDS.id,
        start: this.currentCDS.start,
        stop: this.currentCDS.stop,
        strand: this.currentCDS.strand,
        function: "None",
        status: "Pass"
      };
      this.updateCDS(payload, this.currentCDS.id);
    },
    setFunction(funct) {
      this.newFunction = funct;
    },
    setStart(strt) {
      this.newStart = strt;
    }
  }
};
</script>

<style scoped>
/* ----- Title and Headers ----- */
.headers {
  margin: 40px auto;
}

.alert-primary {
  text-align: left;
  margin: 40px auto;
}

.subheader {
  text-align: left;
}

.info-bottom {
  margin: 50px auto;
}

/* ----- Coding Potential ----- */
.table-responsive {
  max-height: 250px;
  overflow-y: auto;
  display: inline-block;
}

.table-responsive thead th {
  position: sticky;
  top: 0;
  background: #eee;
  border: darkgray;
}

caption {
  display: table-caption;
  /* caption-side: top; */
}

tbody {
  width: 100%;
}

.coding-potential-graphs {
  height: 670px;
}

/* ----- Blast Results ----- */
.btn-blast {
  width: 100%;
}

/* --------------------- */
/* Responsive Design */
@media only screen and (max-width: 1200px) {
  .coding-potential-graphs {
    height: 1300px;
  }
}
</style>
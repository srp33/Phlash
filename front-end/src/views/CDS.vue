<template>
  <div class="container">
    <loading
      :active.sync="pageLoading"
      :is-full-page="true"
      :height="100"
      :width="100"
    ></loading>
    <div class="headers">
      <h1>ID: {{ $route.params.cdsID }}</h1>
      <h4>Start: {{ currentCDS.start }}</h4>
      <h4>Stop: {{ currentCDS.stop }}</h4>
      <h4>Strand: {{ currentCDS.strand }}</h4>
      <h4>Status: {{ currentCDS.status }}</h4>
    </div>
    <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>
        Choose a new start for this gene call based on the information given
        below or keep the current start.
      </p>
      <p>
        <strong>Your selected gene span:</strong> {{ newStart }}-{{ newStop }}
      </p>
      <p><strong>Your selected function:</strong> {{ newFunction }}</p>
      <router-link
        :to="{
          name: 'Annotations',
          params: { phageID: $route.params.phageID },
        }"
      >
        <button class="btn btn-light btn-action">
          <strong>Return</strong>
        </button>
      </router-link>
      <button
        type="button"
        class="btn btn-light btn-action"
        @click="deleteCDS($route.params.cdsID)"
      >
        <strong>Delete</strong>
      </button>
      <button
        type="button"
        v-if="newFunction != ''"
        class="btn btn-light btn-action"
        @click="editCDS"
      >
        <strong>Update</strong>
      </button>
    </div>
    <div class="coding-potential-table">
      <h4 style="text-align: center; margin: 20px">Alternative Gene Spans</h4>
      <div class="table-responsive">
        <table id="cp-table" class="table table-hover">
          <thead>
            <tr>
              <th scope="col">Gene Span</th>
              <th scope="col">Action</th>
            </tr>
          </thead>
          <tbody>
            <tr v-for="(start, index) in startOptions" :key="index">
              <th v-if="start === currentCDS.start">
                {{ start }}-{{ stopOptions[index] }} (current)
              </th>
              <th v-else>{{ start }}-{{ stopOptions[index] }}</th>
              <td>
                <button
                  class="btn btn-dark btn-sm"
                  @click="setGeneSpan(start, stopOptions[index])"
                >
                  <strong>Select</strong>
                </button>
              </td>
            </tr>
          </tbody>
        </table>
      </div>
    </div>
    <hr />
    <div class="coding-potential">
      <h4 style="text-align: center; margin: 40px; height: 100%">
        GeneMark's Coding Potential Per Frame
      </h4>
      <div style="display: inline">
        <div
          style="
            width: 52%;
            display: inline-block;
            float: left;
            margin-right: 10px;
          "
        >
          <strong>Direct Sequences</strong>
        </div>
        <div style="width: 19%; display: inline-block">
          <strong>Complementary Sequences</strong>
        </div>
      </div>
      <div class="coding-potential-graphs">
        <div v-if="dataExists">
          <Graphs
            :data1="data1"
            :data2="data2"
            :data3="data3"
            :data4="data4"
            :data5="data5"
            :data6="data6"
            :start="currentCDS.start"
            :stop="currentCDS.stop"
            :frame="frame"
          />
        </div>
        <div v-else>
          <p>Loading graphs...</p>
        </div>
      </div>
      <p class="graphs-caption">
        <em>Frames are counted 1-6 (direct 1-3 and complementary 4-6).</em>
      </p>
    </div>
    <hr />
    <div class="blast-results">
      <h4 style="text-align: center; margin: 40px">BLAST Results</h4>
      <div v-if="dataExists" id="accordion">
        <div class="card" v-for="key in sortedBlastKeys" :key="key">
          <div class="card-header">
            <h4 class="mb-0">
              <button
                class="btn btn-light btn-blast"
                data-toggle="collapse"
                aria-expanded="false"
                v-bind:data-target="'#' + key"
                v-bind:aria-controls="key"
              >
                <strong v-if="key === currentCDS.start.toString()"
                  >{{ key }} (current)</strong
                >
                <strong v-else>{{ key }}</strong>
              </button>
            </h4>
          </div>
          <div v-bind:id="key" class="collapse" data-parent="#accordion">
            <div class="card-body">
              <BlastResults
                v-if="key === currentCDS.start.toString()"
                :blastResults="blastResults[key]"
                :allowSelect="true"
                @newFunction="setFunction"
              />
              <BlastResults
                v-if="key != currentCDS.start.toString()"
                :blastResults="blastResults[key]"
                :allowSelect="false"
                @newFunction="setFunction"
              />
            </div>
          </div>
        </div>
      </div>
    </div>
    <div class="info-bottom">
      <p>
        <strong>Your selected gene span:</strong> {{ newStart }}-{{ newStop }}
      </p>
      <p><strong>Your function selection:</strong> {{ newFunction }}</p>
      <router-link
        :to="{
          name: 'Annotations',
          params: { phageID: $route.params.phageID },
        }"
      >
        <button class="btn btn-light btn-action">
          <strong>Return</strong>
        </button>
      </router-link>
      <button
        type="button"
        class="btn btn-light btn-action"
        @click="deleteCDS($route.params.cdsID)"
      >
        <strong>Delete</strong>
      </button>
      <button
        type="button"
        v-if="newFunction != ''"
        class="btn btn-light btn-action"
        @click="editCDS"
      >
        <strong>Update</strong>
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
  name: "CDS",
  components: {
    BlastResults,
    Graphs,
    Loading,
  },

  data() {

    return {
      stopOptions: [],
      startOptions: [],
      currentCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: "",
      },
      updatedCDS: {
        id: "",
        start: "",
        stop: "",
        strand: "",
        function: "",
        status: "",
      },
      frame: null,
      blastResults: [],
      showFunction: false,
      showStart: false,
      newFunction: "",
      newStart: null,
      newStop: null,
      dataExists: false,
      pageLoading: true,
      data1: [{ x: [], y: [] }],
      data2: [{ x: [], y: [] }],
      data3: [{ x: [], y: [] }],
      data4: [{ x: [], y: [] }],
      data5: [{ x: [], y: [] }],
      data6: [{ x: [], y: [] }],
      nextCDS: null,
    };

  },

  created() {
    this.getData(this.$route.params.cdsID);
  },

  computed: {

    sortedBlastKeys: function () {
      return Object.keys(this.blastResults).sort((a, b) => b - a);
    },

  },

  methods: {

    getData(cdsID) {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${cdsID}`
        )
        .then((response) => {
          this.currentCDS = response.data.cds;
          this.blastResults = response.data.blast;
          this.startOptions = response.data.start_options;
          this.stopOptions = response.data.stop_options;
          this.newStart = this.currentCDS.start;
          this.newStop = this.currentCDS.stop;
          this.data1 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_1,
            },
          ];
          this.data2 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_2,
            },
          ];
          this.data3 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_3,
            },
          ];
          this.data4 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_4,
            },
          ];
          this.data5 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_5,
            },
          ];
          this.data6 = [
            {
              x: response.data.x_data,
              y: response.data.y_data_6,
            },
          ];
          this.dataExists = true;
          this.pageLoading = false;
          this.frame = ((this.currentCDS.start + 2) % 3) + 1;
          if (this.currentCDS.strand == "-") this.frame += 3;

          this.nextCDS = response.data.nextCDS;
          console.log(this.nextCDS);
        })
        .catch((error) => {
          console.error(error);
        });
    },

    editCDS() {
      this.updatedCDS = this.currentCDS;
      this.updatedCDS.start = this.newStart;
      this.updatedCDS.stop = this.newStop;
      this.updatedCDS.function = this.newFunction;
      const payload = {
        id: this.updatedCDS.id,
        start: this.updatedCDS.start,
        stop: this.updatedCDS.stop,
        strand: this.updatedCDS.strand,
        function: this.updatedCDS.function,
        status: this.currentCDS.status,
      };
      this.updateCDS(payload, this.updatedCDS.id);
    },
    
    updateCDS(payload, cdsID) {
      axios
        .put(
          process.env.VUE_APP_BASE_URL +
            `/annotations/cds/${this.$route.params.phageID}/${cdsID}`,
          payload
        )
        .then(() => {
          if (this.nextCDS != undefined) {
            this.$route.params.cdsID = this.nextCDS;
            this.$router.push(
              `/annotations/cds/${this.$route.params.phageID}/${this.$route.params.cdsID}`
            );
            window.location.reload();
          } else {
            this.$router.push(`/annotations/${this.$route.params.phageID}`);
          }
        })
        .catch((error) => {
          console.error(error);
        });
    },

    deleteCDS(cdsID) {
      this.newFunction = "DELETED";
      this.editCDS();
    },

    keepOriginal() {
      const payload = {
        id: this.currentCDS.id,
        start: this.currentCDS.start,
        stop: this.currentCDS.stop,
        strand: this.currentCDS.strand,
        function: "None",
        status: this.currentCDS.status,
      };
      this.updateCDS(payload, this.currentCDS.id);
    },

    setFunction(funct) {
      this.newFunction = funct;
    },

    setGeneSpan(start, stop) {
      if (start != this.newStart) {
        this.dataExists = false;
        this.newFunction = "";
        this.newStop = stop;
        this.currentCDS.stop = stop;
        this.newStart = start;
        this.currentCDS.start = start;
        this.frame = ((start + 2) % 3) + 1;
        if (this.currentCDS.strand == "-") this.frame += 3;
        this.$nextTick().then(() => {
          this.dataExists = true;
        });
      }
    },

  },
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

.btn-action {
  margin: 7px;
}

/* ----- Coding Potential ----- */
.coding-potential-table {
  margin: 50px auto;
  width: 40%;
}

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
  height: 650px;
}

.graphs-caption {
  margin-top: 15px;
  text-align: left;
  color: grey;
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

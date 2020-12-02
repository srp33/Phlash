<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <h1>Genome Map</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          This is a visual representation of all of the genes in their current state.<br />
          Scroll right to see the rest of the genome.
        </p>
        <p><strong>Key</strong></p>
        <p>
          <strong class="orange-text">Orange:</strong> An arrow pointing left indicates a gene on the complimentary strand.<br />
          <strong class="blue-text">Blue:</strong> An arrow pointing right indicates a gene on the direct strand.<br />
        </p>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-light btn-nav">
              <svg
                class="bi bi-arrow-left"
                width="1em"
                height="1em"
                viewBox="0 0 16 16"
                fill="currentColor"
                xmlns="http://www.w3.org/2000/svg"
              >
                <path
                  fill-rule="evenodd"
                  d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z"
                  clip-rule="evenodd"
                />
                <path
                  fill-rule="evenodd"
                  d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z"
                  clip-rule="evenodd"
                />
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
        </div>
      </div>
    </div>
    <img v-bind:src="image" />
  </div>
</template>

<script>
import axios from "axios";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";
import Navbar from "../components/Navbar.vue";


export default {
  name: "GeneMap",
  components: {
    Loading,
    Navbar,
  },

  data() {

    return {
      pageLoading: true,
      image: null,
    };
  },

  created() {
    this.getGraph();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      return true;
    },

    navGeneMap: function () {
      return true;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },

  },

  methods: {
      getGraph() {
          axios
            .get(
                process.env.VUE_APP_BASE_URL +
                    `/annotations/geneMap/${this.$route.params.phageID}`
            )
            .then((response) => {
                this.image = "data:image/png;base64, " + response.data.image.slice(2, response.data.image.length - 1);
                this.pageLoading = false;
            })
            .catch((error) => {
                console.log(error);
            });
      }
  },
}
</script>

<style scoped>

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 10px;
}

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

.blue-text {
  color: rgb(75, 114, 185);
}

.orange-text {
  color: rgb(233, 139, 105);
}

</style>

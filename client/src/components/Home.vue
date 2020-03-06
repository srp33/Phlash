<template>
<div class="container">

    <h1><strong><i>PHLASH</i></strong></h1>
    <br>
    <div align="left">
        <h2>What is <i>Phlash?</i></h2>
        <p>
            <i>Phlash</i> is an interactive web application that automates the
            manual process of bacteriophage genome annotation by combining
            functionalities of DNA Master and GeneMark and outputting a
            submission-ready GenBank file. By dramatically accelerating the pace
            of bacteriophage genome annotation, <i>Phlash</i> will enable phage
            hunters to spend less time annotating and more time contributing to
            our knowledge on the genetic diversity of bacteriophages and their
            implications. 
        </p>
        <p>
            Please click the button below to get started!
        </p>
        <!-- <ul>
            <li>DNAMaster instructions</li>
            <li>GeneMark instructions</li>
        </ul> -->
    </div>
    <router-link to="/upload">
        <button class="btn btn-primary">Next: Upload DNA Master file</button>
    </router-link>
</div>
</template>

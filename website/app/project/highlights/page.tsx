"use client";

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { motion } from "framer-motion";

export default function HighlightsPage() {
  return (
    <div className="container mx-auto px-4 py-12">
      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6 }}
        className="mb-12"
      >
        <h1 className="mb-4 text-5xl font-bold">
          <span className="bg-gradient-to-r from-primary to-accent bg-clip-text text-transparent">
            Key Highlights
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Selected excerpts from the Ignition Protocol project book.
        </p>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="space-y-8"
      >
        {/* Abstract */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.3 }}
        >
          <Card className="border-primary/50 bg-primary/5">
          <CardHeader>
            <CardTitle className="text-2xl">Abstract</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-muted-foreground">
              Ignition Protocol is a computational framework for discovering, filtering, and operationalising 
              synthetic fuel candidates. In the current implementation, the system is a unified pipeline with two 
              coupled engines:
            </p>
            <ul className="space-y-2 text-muted-foreground">
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span>
                  <strong className="text-foreground">Molecular Generation Module:</strong> Evolves candidate molecules 
                  from fragment libraries using RDKit-based validity checks and a multi-objective fitness score.
                </span>
              </li>
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span>
                  <strong className="text-foreground">Blend Optimization Engine:</strong> Deterministic optimization that 
                  selects and proportions real components into lab-ready fuel blends subject to fuel-type specifications.
                </span>
              </li>
            </ul>
            <p className="text-muted-foreground">
              The project is designed to bridge computational fuel design with proposed experimental validation: 
              combustion "fingerprinting" (flame imaging and spectroscopy) and environmental stability assays 
              (storage stability and biodegradability screening).
            </p>
          </CardContent>
        </Card>
        </motion.div>

        {/* Introduction */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.4 }}
        >
          <Card>
          <CardHeader>
            <CardTitle className="text-2xl">Introduction</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-muted-foreground">
              Synthetic fuels and fuel blends offer a path to reducing reliance on fossil energy while meeting 
              global demand and tightening emissions constraints. Traditional formulation approaches often revolve 
              around modifying known components, which limits exploration of the much larger chemical space of 
              fuel-like molecules.
            </p>
            <p className="text-muted-foreground">
              Ignition Protocol addresses this by combining:
            </p>
            <ul className="space-y-2 text-muted-foreground">
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">Broad exploration:</strong> A fragment-based evolutionary generator that proposes fuel-like molecules</span>
              </li>
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">Operational realism:</strong> A deterministic blend optimiser that respects real fuel constraints and outputs lab-ready mixing recipes</span>
              </li>
            </ul>
            <p className="text-muted-foreground">
              In other words, the project is not just about ranking "interesting molecules"; it is about producing 
              candidate blends that can be mixed, tested, and iterated upon with clear reproducibility.
            </p>
          </CardContent>
        </Card>
        </motion.div>

        {/* Key Methodology Points */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.5 }}
        >
          <Card>
          <CardHeader>
            <CardTitle className="text-2xl">Key Methodology Points</CardTitle>
          </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <h3 className="mb-2 font-semibold text-primary">Fragment-Based Evolution</h3>
              <p className="text-muted-foreground">
                The molecular generator uses fragment libraries (GDB-11 database) to evolve novel fuel molecules 
                through genetic algorithms. The system ensures 100% chemical validity using RDKit validation.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Multi-Objective Fitness</h3>
              <p className="text-muted-foreground">
                Candidates are evaluated through a multi-objective fitness function balancing molecular weight, 
                oxygen-to-carbon ratio, volatility, and combustion potential. The system preserves diversity through 
                Tanimoto similarity-based diversity control.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Deterministic Optimization</h3>
              <p className="text-muted-foreground">
                The blend optimizer uses a deterministic greedy-then-refine algorithm that produces identical results 
                for the same seed, ensuring reproducibility for laboratory work. Optimization typically completes in 
                0.1-0.3 seconds for typical component libraries.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Fuel-Type Specific Constraints</h3>
              <p className="text-muted-foreground">
                The system enforces fuel-type specific constraints including volatility limits (T10/T90, RVP), minimum 
                octane/cetane requirements, density bounds (diesel/jet), and caps on each component's maximum volume fraction.
              </p>
            </div>
          </CardContent>
        </Card>
        </motion.div>

        {/* Notable Results */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.6 }}
        >
          <Card>
            <CardHeader>
              <CardTitle className="text-2xl">Notable Results & Achievements</CardTitle>
            </CardHeader>
          <CardContent className="space-y-4">
            <div>
              <h3 className="mb-2 font-semibold text-primary">Unified Pipeline</h3>
              <p className="text-muted-foreground">
                Successfully integrated molecular generation with blend optimization, creating an end-to-end workflow: 
                Generate → Convert → Optimize → Validate.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Lab-Ready Outputs</h3>
              <p className="text-muted-foreground">
                The system produces volume fractions ready for mixing, with complete lab recipes including mL and grams 
                per 1000mL, enabling immediate experimental validation.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">Reproducibility</h3>
              <p className="text-muted-foreground">
                Deterministic optimization ensures that the same inputs produce identical outputs, with complete parameter 
                capture in RUN.json files for full traceability.
              </p>
            </div>
            <div>
              <h3 className="mb-2 font-semibold text-primary">GDB-11 Integration</h3>
              <p className="text-muted-foreground">
                Integrated support for GDB-11 database with memory-efficient streaming, fuel-relevance filtering (C/H/O only), 
                and preprocessing utilities for large-scale molecular exploration.
              </p>
            </div>
          </CardContent>
        </Card>
        </motion.div>

        {/* Future Directions */}
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.7 }}
        >
          <Card>
            <CardHeader>
              <CardTitle className="text-2xl">Future Directions</CardTitle>
            </CardHeader>
          <CardContent className="space-y-4">
            <p className="text-muted-foreground">
              The project is designed to bridge computational fuel design with experimental validation:
            </p>
            <ul className="space-y-2 text-muted-foreground">
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">Combustion Fingerprinting:</strong> Flame imaging and spectroscopy for experimental validation</span>
              </li>
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">Environmental Stability:</strong> Storage stability and biodegradability screening assays</span>
              </li>
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">ML Integration:</strong> Pathway for hierarchical filtering followed by GPU-accelerated models</span>
              </li>
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">Home Synthesis Mode:</strong> Unique feature for generating molecules optimized for home synthesis feasibility with synthesis route prediction, difficulty ratings (1-10), and yield estimates</span>
              </li>
              <li className="flex items-start gap-2">
                <span className="mt-1 text-primary">•</span>
                <span><strong className="text-foreground">Synthesis Feasibility:</strong> Enhanced synthesis route prediction with reaction types and complexity assessment</span>
              </li>
            </ul>
          </CardContent>
        </Card>
        </motion.div>

        {/* Callout */}
        <motion.div
          initial={{ opacity: 0, scale: 0.95 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.5, delay: 0.8 }}
        >
          <Card className="border-primary/50 bg-gradient-to-r from-ember-dark-red/20 via-background to-ember-dark-red/20">
            <CardContent className="pt-6">
              <p className="text-center text-lg text-muted-foreground">
                <strong className="text-foreground">Full Project Book Available:</strong> The complete Ignition Protocol project book with 
                detailed methodology, literature review, results, and analysis is available at the exhibition stand. 
                Visit us to learn more!
              </p>
            </CardContent>
          </Card>
        </motion.div>
      </motion.div>
    </div>
  );
}


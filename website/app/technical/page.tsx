"use client";

import { motion } from "framer-motion";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Database, Cpu, Zap, Settings, FlaskConical, Code } from "lucide-react";

export default function TechnicalPage() {
  const fuelTypes = [
    {
      type: "Gasoline",
      constraints: [
        "T10: 40-70°C",
        "T90: max 180°C",
        "RVP: max 60 kPa",
        "RON: min 87",
      ],
    },
    {
      type: "Diesel",
      constraints: [
        "Cetane: min 40",
        "Density: 0.82-0.86 g/mL",
        "T90: max 370°C",
      ],
    },
    {
      type: "Jet Fuel",
      constraints: [
        "Density: 0.775-0.84 g/mL",
        "T90: max 300°C",
        "PMI: max 25",
      ],
    },
  ];

  const algorithms = [
    {
      name: "Greedy Forward Selection",
      description: "Iteratively selects best components based on multi-objective scoring",
      details: "Uses weighted combination of RON, cetane, LHV, PMI, and cost",
    },
    {
      name: "Coordinate Descent Refinement",
      description: "Fine-tunes component fractions after initial selection",
      details: "Optimizes volume fractions while respecting constraints and caps",
    },
    {
      name: "Simplex Projection",
      description: "Ensures volume fractions sum to 1.0 and respect component caps",
      details: "Handles delta mode with L1 distance budgets",
    },
    {
      name: "Genetic Algorithm Evolution",
      description: "Evolves novel molecules from fragment libraries",
      details: "Uses crossover, mutation, and selection with diversity preservation",
    },
  ];

  const features = [
    {
      category: "Optimization",
      items: [
        "Multi-objective scoring (RON, cetane, LHV, PMI, cost)",
        "Theme-based optimization (high RON, low PMI, high LHV, low cost)",
        "Delta mode: L1 distance budget optimization around baselines",
        "Novel component integration with configurable caps (max_vol_frac)",
        "Waste credit economics for novel components",
        "Multiple random starts for exploration",
        "Time-budget and evaluation limits",
        "Portfolio generation for diverse solution sets",
      ],
    },
    {
      category: "Molecular Generation",
      items: [
        "Fragment-based evolution from GDB-11 (2.7M+ fragments)",
        "GDB-17 scalability ready (billions of molecules)",
        "RDKit-based chemical validation (100% valid output)",
        "Multi-objective fitness functions (MW, logP, O/C ratio, H-bonding)",
        "Diversity preservation (Tanimoto similarity metrics)",
        "Home synthesis mode: optimized for home synthesis feasibility",
        "Synthesis route prediction with difficulty ratings (1-10)",
        "Oxidation susceptibility prediction (0-10 scale)",
        "Stability testing support (oxidation, phase separation, pH)",
        "Reservoir sampling for uniform random access",
      ],
    },
    {
      category: "Output & Reproducibility",
      items: [
        "Lab-ready recipes (volume fractions, mL, grams per 1000mL)",
        "Complete property analysis",
        "RUN.json with all parameters",
        "Portfolio diversity analysis",
        "CSV exports for all results",
        "Deterministic results (same seed = same output)",
      ],
    },
    {
      category: "Performance",
      items: [
        "Numba JIT compilation for speed",
        "GPU acceleration support (JAX optional)",
        "Memory-efficient streaming for billion-molecule databases",
        "Parallel processing support (joblib)",
        "Optimized for 0.1-0.3s typical runtime",
        "Benchmark profiles: Ultra-snappy (1-2s) to Brutal (45-60min)",
        "Linear scaling with candidate count",
      ],
    },
  ];

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
            Technical Specifications
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Detailed technical capabilities and specifications of the Ignition Protocol system.
        </p>
      </motion.div>

      {/* Fuel Type Specifications */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Supported Fuel Types</h2>
        <div className="grid gap-4 md:grid-cols-3">
          {fuelTypes.map((fuel, index) => (
            <motion.div
              key={fuel.type}
              initial={{ opacity: 0, y: 30 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.5, delay: 0.3 + index * 0.1 }}
            >
              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center gap-2">
                    <FlaskConical className="h-5 w-5 text-primary" />
                    {fuel.type}
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <ul className="space-y-2">
                    {fuel.constraints.map((constraint, idx) => (
                      <li key={idx} className="flex items-start gap-2 text-sm text-muted-foreground">
                        <span className="mt-1 text-primary">•</span>
                        <span>{constraint}</span>
                      </li>
                    ))}
                  </ul>
                </CardContent>
              </Card>
            </motion.div>
          ))}
        </div>
      </motion.section>

      {/* Algorithms */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.6 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Core Algorithms</h2>
        <div className="grid gap-4 md:grid-cols-2">
          {algorithms.map((algo, index) => (
            <motion.div
              key={algo.name}
              initial={{ opacity: 0, x: -20 }}
              animate={{ opacity: 1, x: 0 }}
              transition={{ duration: 0.5, delay: 0.7 + index * 0.1 }}
            >
              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center gap-2">
                    <Code className="h-5 w-5 text-primary" />
                    {algo.name}
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <p className="mb-2 text-muted-foreground">{algo.description}</p>
                  <p className="text-sm text-muted-foreground/80">{algo.details}</p>
                </CardContent>
              </Card>
            </motion.div>
          ))}
        </div>
      </motion.section>

      {/* Features by Category */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 1.0 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Feature Overview</h2>
        <div className="grid gap-6 md:grid-cols-2">
          {features.map((category, index) => (
            <motion.div
              key={category.category}
              initial={{ opacity: 0, y: 30 }}
              animate={{ opacity: 1, y: 0 }}
              transition={{ duration: 0.5, delay: 1.1 + index * 0.1 }}
            >
              <Card>
                <CardHeader>
                  <CardTitle>{category.category}</CardTitle>
                </CardHeader>
                <CardContent>
                  <ul className="space-y-2">
                    {category.items.map((item, idx) => (
                      <li key={idx} className="flex items-start gap-2 text-sm text-muted-foreground">
                        <span className="mt-1 text-primary">•</span>
                        <span>{item}</span>
                      </li>
                    ))}
                  </ul>
                </CardContent>
              </Card>
            </motion.div>
          ))}
        </div>
      </motion.section>

      {/* Advanced Features */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 1.2 }}
        className="mb-12"
      >
        <h2 className="mb-6 text-3xl font-bold">Advanced Features</h2>
        <div className="grid gap-4 md:grid-cols-2">
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <FlaskConical className="h-5 w-5 text-primary" />
                Home Synthesis Mode
              </CardTitle>
            </CardHeader>
            <CardContent>
              <p className="mb-3 text-muted-foreground">
                Unique feature for generating molecules optimized for home synthesis feasibility. Includes:
              </p>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• Synthesis feasibility scoring (0-100)</li>
                <li>• Complexity rating (1-10 scale)</li>
                <li>• Suggested synthesis routes with reaction types</li>
                <li>• Yield estimates and difficulty ratings</li>
                <li>• Home-synthesizable fragment library</li>
              </ul>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Database className="h-5 w-5 text-primary" />
                Database Scalability
              </CardTitle>
            </CardHeader>
            <CardContent>
              <p className="mb-3 text-muted-foreground">
                Designed to scale from demo-sized datasets to billion-molecule databases:
              </p>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• Current: GDB-11 (2.7M+ fragments, 11 atoms max)</li>
                <li>• Ready: GDB-17 (billions of molecules, 17 atoms max)</li>
                <li>• Memory-efficient streaming architecture</li>
                <li>• Reservoir sampling for uniform access</li>
                <li>• Preprocessing utilities for fuel-relevance filtering</li>
              </ul>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Zap className="h-5 w-5 text-primary" />
                Experimental Validation Framework
              </CardTitle>
            </CardHeader>
            <CardContent>
              <p className="mb-3 text-muted-foreground">
                Planned experimental protocols for validating generated candidates:
              </p>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• Combustion fingerprinting (flame imaging, spectroscopy)</li>
                <li>• Biodegradability screening (soil slurry assays)</li>
                <li>• Storage stability tests (oxidation, phase separation)</li>
                <li>• Feature extraction and PCA visualization</li>
              </ul>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Settings className="h-5 w-5 text-primary" />
                Benchmark Profiles
              </CardTitle>
            </CardHeader>
            <CardContent>
              <p className="mb-3 text-muted-foreground">
                Scalable performance profiles from instant feedback to deep exploration:
              </p>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• Ultra-snappy: 30-40 candidates, 1-2 seconds</li>
                <li>• Demo: 50-100 candidates, 5-10 seconds</li>
                <li>• Classroom: 100-200 candidates, 20-50 generations</li>
                <li>• Brutal: 800,000+ candidates, 45-60 minutes</li>
              </ul>
            </CardContent>
          </Card>
        </div>
      </motion.section>

      {/* Technical Stack */}
      <motion.section
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 1.4 }}
      >
        <h2 className="mb-6 text-3xl font-bold">Technical Stack</h2>
        <div className="grid gap-4 md:grid-cols-3">
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Database className="h-5 w-5 text-primary" />
                Core Libraries
              </CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• NumPy, Pandas, SciPy</li>
                <li>• RDKit (cheminformatics)</li>
                <li>• Numba (JIT compilation)</li>
                <li>• scikit-optimize</li>
              </ul>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Zap className="h-5 w-5 text-primary" />
                Performance
              </CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• JAX (GPU acceleration)</li>
                <li>• Joblib (parallel processing)</li>
                <li>• Memory-efficient streaming</li>
                <li>• Optimized algorithms</li>
              </ul>
            </CardContent>
          </Card>
          <Card>
            <CardHeader>
              <CardTitle className="flex items-center gap-2">
                <Settings className="h-5 w-5 text-primary" />
                Data Formats
              </CardTitle>
            </CardHeader>
            <CardContent>
              <ul className="space-y-1 text-sm text-muted-foreground">
                <li>• CSV (components, blends)</li>
                <li>• JSON (metadata, config)</li>
                <li>• SMILES (molecules)</li>
                <li>• GDB-11 (fragments)</li>
              </ul>
            </CardContent>
          </Card>
        </div>
      </motion.section>
    </div>
  );
}


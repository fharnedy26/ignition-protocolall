"use client";

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { motion } from "framer-motion";

export default function ExamplesPage() {
  const examples = [
    {
      title: "Basic Gasoline Optimization",
      description: "Simple optimization with default parameters",
      code: `python engine.py \\
  --components data/components.csv \\
  --fuel-type gasoline \\
  --K 5 \\
  --top 10`,
    },
    {
      title: "Delta Mode (Small Changes)",
      description: "Make small edits around an existing blend",
      code: `python engine.py \\
  --components data/components.csv \\
  --baseline data/baseline_blend.csv \\
  --fuel-type gasoline \\
  --K 5 \\
  --top 15 \\
  --delta-budget 0.2`,
    },
    {
      title: "With Novel Components",
      description: "Include novel components with a cap",
      code: `python engine.py \\
  --components data/components.csv \\
  --novel data/components_novel.csv \\
  --allow-novel 1 \\
  --max-add-novel 0.1 \\
  --fuel-type gasoline \\
  --K 5 \\
  --top 10`,
    },
    {
      title: "Diesel Optimization",
      description: "Optimize for diesel fuel specifications",
      code: `python engine.py \\
  --components data/components.csv \\
  --fuel-type diesel \\
  --K 4 \\
  --top 15 \\
  --seed 123`,
    },
    {
      title: "Jet Fuel Optimization",
      description: "Optimize for jet fuel specifications",
      code: `python engine.py \\
  --components data/components.csv \\
  --fuel-type jet \\
  --K 5 \\
  --top 20 \\
  --seed 456`,
    },
    {
      title: "Using Presets",
      description: "Use predefined configuration presets",
      code: `# Baseline preset (strict specs, no novel)
python engine.py --components data/components.csv \\
  --preset baseline --fuel-type gasoline --top 5

# Delta preset (small changes around baseline)
python engine.py --components data/components.csv \\
  --baseline data/baseline_blend.csv \\
  --preset delta --fuel-type gasoline --top 5

# Novel preset (allow up to 3% novel components)
python engine.py --components data/components.csv \\
  --novel data/components_novel.csv \\
  --preset novel3 --fuel-type gasoline --top 5`,
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
            Usage Examples
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Real-world usage scenarios and code samples for different optimization tasks.
        </p>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="space-y-8"
      >
        {examples.map((example, idx) => (
          <motion.div
            key={idx}
            initial={{ opacity: 0, y: 30 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5, delay: 0.3 + idx * 0.1 }}
          >
            <Card>
            <CardHeader>
              <CardTitle>{example.title}</CardTitle>
              <p className="text-muted-foreground">{example.description}</p>
            </CardHeader>
            <CardContent>
              <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                <code>{example.code}</code>
              </pre>
            </CardContent>
          </Card>
          </motion.div>
        ))}

        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.8 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Python API Examples</h2>
          <Card>
            <CardHeader>
              <CardTitle>Programmatic Usage</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
              <div>
                <h3 className="mb-2 font-semibold">Basic Optimization</h3>
                <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                  <code>{`from fuel_engine import greedy_then_refine, load_components

comps = load_components("data/components.csv")
results = greedy_then_refine(comps, K=5, fuel_type="gasoline")

# Get best result
best = results[0]
print(f"Best score: {best['score']}")
for comp in best['components']:
    print(f"{comp['id']}: {comp['vol_frac']:.3f}")`}</code>
                </pre>
              </div>
              <div>
                <h3 className="mb-2 font-semibold">With Baseline</h3>
                <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                  <code>{`from fuel_engine import greedy_then_refine, load_components, load_baseline

comps = load_components("data/components.csv")
baseline = load_baseline("data/baseline_blend.csv")

results = greedy_then_refine(
    comps,
    baseline=baseline,
    K=5,
    fuel_type="gasoline",
    delta_budget=0.2
)`}</code>
                </pre>
              </div>
            </CardContent>
          </Card>
        </motion.section>
      </motion.div>
    </div>
  );
}


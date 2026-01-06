"use client";

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { motion } from "framer-motion";

export default function APIPage() {
  const functions = [
    {
      name: "greedy_then_refine",
      description: "Main optimization function that performs greedy selection followed by refinement",
      params: [
        { name: "comps", type: "DataFrame", desc: "Component library DataFrame" },
        { name: "K", type: "int", desc: "Maximum components in blend" },
        { name: "fuel_type", type: "str", desc: "Fuel type: 'gasoline', 'diesel', or 'jet'" },
        { name: "strict", type: "bool", desc: "Enforce strict specifications" },
        { name: "baseline", type: "dict", desc: "Optional baseline blend dictionary" },
        { name: "allow_novel", type: "bool", desc: "Allow novel components" },
        { name: "max_novel", type: "float", desc: "Maximum novel fraction (0-1)" },
      ],
      returns: "List[Dict] - List of optimized blend results",
    },
    {
      name: "load_components",
      description: "Load component library from CSV file",
      params: [
        { name: "file_path", type: "str", desc: "Path to component CSV file" },
      ],
      returns: "DataFrame - Component library with validated columns",
    },
    {
      name: "blend_props",
      description: "Calculate blend properties from component fractions",
      params: [
        { name: "x", type: "ndarray", desc: "Volume fractions array" },
        { name: "comps", type: "DataFrame", desc: "Component library" },
      ],
      returns: "Dict - Blend properties (RON, MON, cetane, LHV, PMI, TSI, etc.)",
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
            API Reference
          </span>
        </h1>
        <p className="text-xl text-muted-foreground">
          Complete API documentation for the fuel blend optimization engine.
        </p>
      </motion.div>

      <motion.div
        initial={{ opacity: 0, y: 20 }}
        animate={{ opacity: 1, y: 0 }}
        transition={{ duration: 0.6, delay: 0.2 }}
        className="space-y-8"
      >
        {functions.map((func, index) => (
          <motion.div
            key={func.name}
            initial={{ opacity: 0, y: 30 }}
            animate={{ opacity: 1, y: 0 }}
            transition={{ duration: 0.5, delay: 0.3 + index * 0.1 }}
          >
            <Card>
            <CardHeader>
              <CardTitle className="font-mono text-2xl">{func.name}</CardTitle>
              <p className="text-muted-foreground">{func.description}</p>
            </CardHeader>
            <CardContent className="space-y-4">
              <div>
                <h3 className="mb-2 font-semibold">Parameters</h3>
                <div className="space-y-2">
                  {func.params.map((param) => (
                    <div key={param.name} className="rounded-lg bg-muted p-3">
                      <div className="flex items-start gap-2">
                        <code className="font-mono text-primary">{param.name}</code>
                        <span className="text-muted-foreground">:</span>
                        <code className="font-mono text-accent">{param.type}</code>
                      </div>
                      <p className="mt-1 text-sm text-muted-foreground">{param.desc}</p>
                    </div>
                  ))}
                </div>
              </div>
              <div>
                <h3 className="mb-2 font-semibold">Returns</h3>
                <code className="rounded-lg bg-muted p-3 font-mono text-sm">{func.returns}</code>
              </div>
            </CardContent>
          </Card>
          </motion.div>
        ))}

        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.6 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Module Structure</h2>
          <Card>
            <CardContent className="pt-6">
              <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                <code>{`from fuel_engine import (
    greedy_then_refine,
    load_components,
    load_baseline,
    blend_props,
    write_top_n_csv,
    write_run_json,
)
from fuel_engine.constants import SPECS`}</code>
              </pre>
            </CardContent>
          </Card>
        </motion.section>

        <motion.section
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5, delay: 0.8 }}
        >
          <h2 className="mb-4 text-3xl font-bold">Example Usage</h2>
          <Card>
            <CardContent className="pt-6">
              <pre className="overflow-x-auto rounded-lg bg-muted p-4 text-sm">
                <code>{`from fuel_engine import greedy_then_refine, load_components

# Load components
comps = load_components("data/components.csv")

# Optimize
results = greedy_then_refine(
    comps,
    K=5,
    fuel_type="gasoline",
    strict=True
)

# Access best result
best = results[0]
print(f"Score: {best['score']}")
print(f"Components: {best['components']}")`}</code>
              </pre>
            </CardContent>
          </Card>
        </motion.section>
      </motion.div>
    </div>
  );
}

